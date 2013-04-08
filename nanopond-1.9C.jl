# *********************************************************************** #
#                                                                         #
# Nanopond version 1.9 -- A teeny tiny artificial life virtual machine    #
# Copyright (C) 2005 Adam Ierymenko - http://www.greythumb.com/people/api #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program; if not, write to the Free Software             #
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110 USA     #
#                                                                         #
# *********************************************************************** #

#
# Changelog:
#
# 1.0 - Initial release
# 1.1 - Made empty cells get initialized with 0xffff... instead of zeros
#       when the simulation starts. This makes things more consistent with
#       the way the output buf is treated for self-replication, though
#       the initial state rapidly becomes irrelevant as the simulation
#       gets going.  Also made one or two very minor performance fixes.
# 1.2 - Added statistics for execution frequency and metabolism, and made
#       the visualization use 16bpp color.
# 1.3 - Added a few other statistics.
# 1.4 - Replaced SET with KILL and changed EAT to SHARE. The SHARE idea
#       was contributed by Christoph Groth (http://www.falma.de/). KILL
#       is a variation on the original EAT that is easier for cells to
#       make use of.
# 1.5 - Made some other instruction changes such as XCHG and added a
#       penalty for failed KILL attempts. Also made access permissions
#       stochastic.
# 1.6 - Made cells all start facing in direction 0. This removes a bit
#       of artificiality and requires cells to evolve the ability to
#       turn in various directions in order to reproduce in anything but
#       a straight line. It also makes pretty graphics.
# 1.7 - Added more statistics, such as original lineage, and made the
#       genome dump files CSV files as well.
# 1.8 - Fixed LOOP/REP bug reported by user Sotek.  Thanks!  Also
#       reduced the default mutation rate a bit.
# 1.9 - Added a bunch of changes suggested by Christoph Groth: a better
#       coloring algorithm, right click to switch coloring schemes (two
#       are currently supported), and a few speed optimizations. Also
#       changed visualization so that cells with generations less than 2
#       are no longer shown.

# 1.9B - Added a bunch of changes and features.
#
# Changes: (1) Potential speed increases by using bzero(), memset(),
# and memcpy() instead of loops. (Is this really faster? I haven't
# tested it.) (2) added instructions below for compiling on Mac
# OSX. (3) Some code reshuffling, e.g. moved cell execution loop to
# executeCell() function. (4) Removed "facing" variable; cell now
# gets direction from register. (5) Moved main loop into a doMainLoop()
# function.
#
# These new "features" can be turned off by unsetting USE_EXTRAS: (1)
# New keyboard inputs and display settings, e.g. =/- increase/
# decrease avg. energy input. (2) Ability to save/read current run
# to/from a file, and do screenshots to a BMP file. (3) Inserting a
# genome into a cell (by typing in a string) (4) New display modes:
# LOGO and ENERGY (ENERGY is not working correctly - I want it to
# eventually be a greyscale). (5) Replaced TURN command with EJECT,
# eject output buffer into neighbor, potentially at a non-zero genome
# location (could allow for parasite/viral/sexual types of
# behavior). (6) Can add/remove "barriers" to divide the landscape
# (contintental drift!)  (7) 4 STOPs in a row are required for ending
# execution (changed from one; thus it can also act as a NOP). (8)
# In-place mutations of individual living cells, in addition to
# replacement with random genomes.
#/

#
# Nanopond is just what it says: a very very small and simple artificial
# life virtual machine.
#
# It is a "small evolving program" based artificial life system of the same
# general class as Tierra, Avida, and Archis.  It is written in very tight
# and efficient C code to make it as fast as possible, and is so small that
# it consists of only one .c file.
#
# How Nanopond works:
#
# The Nanopond world is called a "pond."  It is an NxN two dimensional
# array of Cell structures, and it wraps at the edges (it's toroidal).
# Each Cell structure consists of a few attributes that are there for
# statistics purposes, an energy level, and an array of POND_DEPTH
# four-bit values.  (The four-bit values are actually stored in an array
# of machine-size words.)  The array in each cell contains the genome
# associated with that cell, and POND_DEPTH is therefore the maximum
# allowable size for a cell genome.
#
# The first four bit value in the genome is called the "logo." What that is
# for will be explained later. The remaining four bit values each code for
# one of 16 instructions. Instruction zero (0x0) is NOP (no operation) and
# instruction 15 (0xf) is STOP (stop cell execution). Read the code to see
# what the others are. The instructions are exceptionless and lack fragile
# operands. This means that *any* arbitrary sequence of instructions will
# always run and will always do *something*. This is called an evolvable
# instruction set, because programs coded in an instruction set with these
# basic characteristics can mutate. The instruction set is also
# Turing-complete, which means that it can theoretically do anything any
# computer can do. If you're curious, the instruciton set is based on this:
# http://www.muppetlabs.com/~breadbox/bf/  (Brainfuck)
#
# At the center of Nanopond is a core loop. Each time this loop executes,
# a clock counter is incremented and one or more things happen:
#
# - Every REPORT_FREQUENCY clock ticks a line of comma seperated output
#   is printed to STDOUT with some statistics about what's going on.
# - Every DUMP_FREQUENCY clock ticks, all viable replicators (cells whose
#   generation is >= 2) are dumped to a file on disk.
# - Every INFLOW_FREQUENCY clock ticks a random x,y location is picked,
#   energy is added (see INFLOW_RATE_MEAN and INFLOW_RATE_DEVIATION)
#   and it's genome is filled with completely random bits.  Statistics
#   are also reset to generation==0 and parentID==0 and a new cell ID
#   is assigned.
# - Every tick a random x,y location is picked and the genome inside is
#   executed until 4 STOP instructions are encountered or the cell's
#   energy counter reaches zero. (Each instruction costs one unit energy.)
#
# The cell virtual machine is an extremely simple register machine with
# a single four bit register, one memory pointer, one spare memory pointer
# that can be exchanged with the main one, and an output buffer. When
# cell execution starts, this output buffer is filled with all binary 1's
# (0xffff....). When cell execution is finished, if the first byte of
# this buffer is *not* 0xff, then the VM says "hey, it must have some
# data!". This data is a candidate offspring; to reproduce cells must
# copy their genome data into the output buffer.
#
# When the VM sees data in the output buffer, it looks at the cell
# adjacent to the cell that just executed and checks whether or not
# the cell has permission (see below) to modify it. If so, then the
# contents of the output buffer replace the genome data in the
# adjacent cell. Statistics are also updated: parentID is set to the
# ID of the cell that generated the output and generation is set to
# one plus the generation of the parent.
#
# A cell is permitted to access a neighboring cell if:
#    - That cell's energy is zero
#    - That cell's parentID is zero
#    - That cell's logo (remember?) matches the trying cell's "guess"
#
# Since randomly introduced cells have a parentID of zero, this allows
# real living cells to always replace them or eat them.
#
# The "guess" is merely the value of the register at the time that the
# access attempt occurs.
#
# Permissions determine whether or not an offspring can take the place
# of the contents of a cell and also whether or not the cell is allowed
# to EAT (an instruction) the energy in it's neighbor.
#
# If you haven't realized it yet, this is why the final permission
# criteria is comparison against what is called a "guess." In conjunction
# with the ability to "eat" neighbors' energy, guess what this permits?
#
# Since this is an evolving system, there have to be mutations. The
# MUTATION_RATE sets their probability. Mutations are random variations
# with a frequency defined by the mutation rate to the state of the
# virtual machine while cell genomes are executing. Since cells have
# to actually make copies of themselves to replicate, this means that
# these copies can vary if mutations have occurred to the state of the
# VM while copying was in progress.
#
# What results from this simple set of rules is an evolutionary game of
# "corewar." In the beginning, the process of randomly generating cells
# will cause self-replicating viable cells to spontaneously emerge. This
# is something I call "random genesis," and happens when some of the
# random gak turns out to be a program able to copy itself. After this,
# evolution by natural selection takes over. Since natural selection is
# most certainly *not* random, things will start to get more and more
# ordered and complex (in the functional sense). There are two commodities
# that are scarce in the pond: space in the NxN grid and energy. Evolving
# cells compete for access to both.
#
# If you want more implementation details such as the actual instruction
# set, read the source. It's well commented and is not that hard to
# read. Most of it's complexity comes from the fact that four-bit values
# are packed into machine size words by bit shifting. Once you get that,
# the rest is pretty simple.
#
# Nanopond, for it's simplicity, manifests some really interesting
# evolutionary dynamics. While I haven't run the kind of multiple-
# month-long experiment necessary to really see this (I might!), it
# would appear that evolution in the pond doesn't get "stuck" on just
# one or a few forms the way some other simulators are apt to do.
# I think simplicity is partly reponsible for this along with what
# biologists call embeddedness, which means that the cells are a part
# of their own world.
#
# Run it for a while... the results can be... interesting!
#
# Running Nanopond:
#
# Nanopond can use SDL (Simple Directmedia Layer) for screen output. If
# you don't have SDL, comment out USE_SDL below and you'll just see text
# statistics and get genome data dumps. (Turning off SDL will also speed
# things up slightly.)
#
# After looking over the tunable parameters below, compile Nanopond and
# run it. Here are some example compilation commands from Linux:
#
# For Pentiums:
#  gcc -O6 -march=pentium -funroll-loops -fomit-frame-pointer -s
#   -o nanopond nanopond.c -lSDL
#
# For Athlons with gcc 4.0+:
#  gcc -O6 -msse -mmmx -march=athlon -mtune=athlon -ftree-vectorize
#   -funroll-loops -fomit-frame-pointer -o nanopond nanopond.c -lSDL
#
# The second line is for gcc 4.0 or newer and makes use of GCC's new
# tree vectorizing feature. This will speed things up a bit by
# compiling a few of the loops into MMX/SSE instructions.
#
# On Mac OS X, I have done it using Fink's version of SDL:
#  gcc -fast -floop-optimize -finline-functions -fstrength-reduce
#     -o nanopond nanopond-1.9.c -D_MSC_VER
#     `/sw/bin/sdl-config --prefix=/sw --libs --cflags`
#
# This should also work on other Posix-compliant OSes with relatively
# new C compilers. (Really old C compilers will probably not work.)
# On other platforms, you're on your own! On Windows, you will probably
# need to find and download SDL if you want pretty graphics and you
# will need a compiler. MinGW and Borland's BCC32 are both free. I
# would actually expect those to work better than Microsoft's compilers,
# since MS tends to ignore C/C++ standards. If stdint.h isn't around,
# you can fudge it like this:
#
# #define uintptr_t unsigned long (or whatever your machine size word is)
# #define uint8_t unsigned char
# #define uint16_t unsigned short
# #define uint64_t unsigned long long (or whatever is your 64-bit int)
#
# When Nanopond runs, comma-seperated stats (see doReport() for
# the columns) are output to stdout and various messages are output
# to stderr. For example, you might do:
#
# ./nanopond >>stats.csv 2>messages.txt &
#
# To get both in seperate files.
#
# <plug>
# Be sure to visit http://www.greythumb.com/blog for your dose of
# technobiology related news. Also, if you're ever in the Boston
# area, visit http://www.greythumb.com/bostonalife to find out when
# our next meeting is!
# </plug>
#
# Have fun!
#/

stop() = error( "STOPPING ON PURPOSE" )

# ----------------------------------------------------------------------- #
# Tunable parameters                                                      #
# ----------------------------------------------------------------------- #

# Iteration to stop at. Comment this out to run forever. #
# global const STOP_AT 150000000000ULL #

# Size of pond in X and Y dimensions. #
global const POND_SIZE_X = 256 ##640
global const POND_SIZE_Y = 256 ##480

# Frequency of comprehensive reports-- lower values will provide more
# info while slowing down the simulation. Higher values will give less
# frequent updates. #
global const REPORT_FREQUENCY = ( 20 * POND_SIZE_X * POND_SIZE_Y ) #10000000#

# Frequency at which to dump all viable replicators (generation > 2)
# to a file named <clock>.dump in the current directory.  Comment
# out to disable. The cells are dumped in hexadecimal, which is
# semi-human-readable if you look at the big switch() statement
# in the main loop to see what instruction is signified by each
# four-bit value. #
global const DUMP_FREQUENCY = ( 10 * REPORT_FREQUENCY ) #100000000#

# Mutation rate -- range is from 0 (none) to 0xffffffff (all mutations!) #
# To get it from a float probability from 0.0 to 1.0, multiply it by
# 4294967295 (0xffffffff) and round. #
global const MUTATION_RATE = 21475 # p=~0.000005 #

# How frequently should random cells / energy be introduced?
# Making this too high makes things very chaotic. Making it too low
# might not introduce enough energy. #
global const INFLOW_FREQUENCY = 100

# Base amount of energy to introduce per INFLOW_FREQUENCY ticks #
global const INFLOW_RATE_BASE = 4000 # About enough to go through 8 genome lengths #

# A random amount of energy between 0 and this is added to
# INFLOW_RATE_BASE when energy is introduced. Comment this out for
# no variation in inflow rate. #
global const INFLOW_RATE_VARIATION = 8000

# Depth of pond in four-bit codons -- this is the maximum
# genome size. This *must* be a multiple of 16! #
global POND_DEPTH = 512

# This is the divisor that determines how much energy is taken
# from cells when they try to KILL a viable cell neighbor and
# fail. Higher numbers mean lower penalties. #
global const FAILED_KILL_PENALTY = 2

# Define this to use SDL. To use SDL, you must have SDL headers
# available and you must link with the SDL library when you compile. #
# Comment this out to compile without SDL visualization support. #
global const USE_SDL = 0

# This is the frequency of screen refreshes if SDL is enabled. #
global const SDL_UPDATE_FREQUENCY = (POND_SIZE_X * POND_SIZE_Y / 2) #200000#

# Use this to turn off all the "extras" that I added and make code very
#   compact again. #
global const USE_EXTRAS = 1

#ifdef USE_EXTRAS
# Maximum allowed number of barriers.
# Undefine this to turn off barriers (and speed it up a bit). #
global const NUM_BARRIERS = 5

# This is the divisor that determines how much energy is taken
# from cells when they try to EJECT into a viable cell neighbor and
# fail. Higher numbers mean lower penalties. #
global const FAILED_EJECT_PENALTY = 2

global const METABOLITE_AMOUNT = 25

# MULTI THREADING! #
#undef N_THREADS
#global const N_THREADS 1

# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
# This is the Mersenne Twister by Makoto Matsumoto and Takuji Nishimura   #
# http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html  #
# ----------------------------------------------------------------------- #

# A few very minor changes were made by me - Adam #

# 
#    A C-program for MT19937, with initialization improved 2002/1/26.
#    Coded by Takuji Nishimura and Makoto Matsumoto.

#    Before using, initialize the state by using init_genrand(seed)  
#    or init_by_array(init_key, key_length).

#    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
#    All rights reserved.                          

#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:

#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.

#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in the
#         documentation and/or other materials provided with the distribution.

#      3. The names of its contributors may not be used to endorse or promote 
#         products derived from this software without specific prior written 
#         permission.

#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
#    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#    Any feedback is very welcome.
#    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
#    email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
#

# Period parameters #  
global const N = 624
global const M = 397
global const MATRIX_A = 0x9908b0df   # constant vector a #
global const UPPER_MASK = 0x80000000 # most significant w-r bits #
global const LOWER_MASK = 0x7fffffff # least significant r bits #

global mt = Array(Uint64, N); # the array for the state vector  #
global mti = N+1; # mti==N+1 means mt[N] is not initialized #

# initializes mt[N] with a seed #
function init_genrand(s::Uint64)
    mt[0]= s & 0xffffffffUL;
    for mti=2:N
        mt[mti] = (1812433253 * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        # See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. #
        # In the previous versions, MSBs of the seed affect   #
        # only MSBs of the array mt[].                        #
        # 2002/01/09 modified by Makoto Matsumoto             #
        mt[mti] &= 0xffffffff;
        # for >32 bit machines #
    end
end

# generates a random number on [0,0xffffffff]-interval #
function genrand_int32()
    y::Uint32;
    mag01::Array{Uint32,1} = [0x0, MATRIX_A];
    # mag01[x] = x * MATRIX_A  for x=0,1 #

    if mti >= N # generate N words at one time #
        for kk=1:N-M
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        end
        for kk=N-M:N
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        end
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    end
  
    y = mt[mti++];

    # Tempering #
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    y;
end

# ----------------------------------------------------------------------- #

# Pond depth in machine-size words.  This is calculated from
# POND_DEPTH and the size of the machine word. (The multiplication
# by two is due to the fact that there are two four-bit values in
# each eight-bit byte.) #
#global const POND_DEPTH_SYSWORDS = (POND_DEPTH / (sizeof(uintptr_t) * 2))
#global const POND_DEPTH_BYTES (POND_DEPTH / 2)

# Number of bits in a machine-size word #
#global const SYSWORD_BITS (sizeof(uintptr_t) * 8)

# Constants representing neighbors in the 2D grid. #
global const N_LEFT = 0
global const N_RIGHT = 1
global const N_UP = 2
global const N_DOWN = 3

# Word and bit at which to start execution #
# This is after the "logo" #
#global const EXEC_START_WORD 0
global const EXEC_START_IND = 2 #4

# Number of bits set in binary numbers 0000 through 1111 #
global const BITS_IN_FOURBIT_WORD::Vector{Uint64} = uint64([ 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 ]);

#*
# Structure for a cell in the pond
#/
type Cell
    # Globally unique cell ID #
    ID::Uint64;
    
    # ID of the cell's parent #
    parentID::Uint64;
    
    # Counter for original lineages -- equal to the cell ID of
    # the first cell in the line. #
    lineage::Uint64;
    
    # Generations start at 0 and are incremented from there. #
    generation::Uint64;

    # Count number of children this cell has spawned. #
    children::Uint64;
    
    # Energy level of this cell #
    energy::Uint64; ##internalEnergy::Uint64;

    # Memory space for cell genome (genome is stored as four
    # bit instructions packed into machine size words) #
    genome::Vector{Uint8}; ##[POND_DEPTH_SYSWORDS];
end;

Cell() = Cell(uint64(0), uint64(0), uint64(0), uint64(0), uint64(0), uint64(0), zeros(Uint8,POND_DEPTH));

# The pond is a 2D array of cells #
##struct Cell pond[POND_SIZE_X][POND_SIZE_Y];
pond = [ Cell() for i=1:POND_SIZE_X, j=1:POND_SIZE_Y ];

#*
# Structure for VM buffers and stacks
#/
type VMStruct
  # Buffer used for execution output of candidate offspring #
  outputBuf::Vector{Uint64}; ##[POND_DEPTH_SYSWORDS];

  # Virtual machine loop/rep stack pointer #
  loopStackPtr::Uint64;
   
  # Virtual machine loop/rep stack #   
  loopStack_wordPtr::Vector{Uint64}; ##[POND_DEPTH];
  loopStack_shiftPtr::Vector{Uint64}; ##[POND_DEPTH];
end

VMStruct() = VMStruct(zeros(Uint64, POND_DEPTH), uint64(0), zeros(Uint64, POND_DEPTH), 
                      zeros(Uint64, POND_DEPTH));

# Currently selected color scheme #
KINSHIP = 1;
LOGO = 2;
LINEAGE = 3;
ENERGY = 4;
MAX_COLOR_SCHEME = 5;
colorScheme = KINSHIP;

##enum { KINSHIP,LOGO,LINEAGE,ENERGY,MAX_COLOR_SCHEME } colorScheme = KINSHIP;
colorSchemeName = { "KINSHIP", "LOGO", "LINEAGE", "ENERGY" };
colorOnlyViable = 0;

#ifdef USE_EXTRAS
instructionNames = { "ZERO", "FWD", "BACK", "INC", "DEC", "READG", "WRITEG", "READB",
		    "WRITEB", "LOOP", "REP", "XCHG", "EJECT", "KILL", "SHARE", "STOP" };
tracePtr = Nothing; ## ::Cell
#ifdef NUM_BARRIERS
static uintptr_t barriers_x[ NUM_BARRIERS ], barriers_y[ NUM_BARRIERS ], n_barriers = 0;
#endif
#endif

# Clock is incremented on each core loop #
static uint64_t clockTicks = 0;
  
# This is used to generate unique cell IDs #
static uint64_t cellIdCounter = 0;

# This is used allow the user to increase or decrease mean energy input #
static intptr_t adjustedMeanEnergy = 0;

static uintptr_t paused = 0;

#*
# Get a random number
#
# @return Random number
#/
static inline uintptr_t getRandom()
{
  # A good optimizing compiler should optimize out this if #
  # This is to make it work on 64-bit boxes #
  if (sizeof(uintptr_t) == 8)
    return (uintptr_t)((((uint64_t)genrand_int32()) << 32) ^ ((uint64_t)genrand_int32()));
  else return (uintptr_t)genrand_int32();
}

#*
# Structure for keeping some running tally type statistics
#/
struct PerReportStatCounters
{
  # Counts for the number of times each instruction was
   * executed since the last report. #
  double instructionExecutions[16];
  
  # Number of cells executed since last report #
  double cellExecutions;
  
  # Number of viable cells replaced by other cells' offspring #
  uintptr_t viableCellsReplaced;
  
  # Number of viable cells KILLed #
  uintptr_t viableCellsKilled;
  
  # Number of successful SHARE operations #
  uintptr_t viableCellShares;

  # Number of viable cells that had code injected into them at a non-zero location #
  uintptr_t bababooms;
};

# Global statistics counters #
struct PerReportStatCounters statCounters;

#*
# Function to decide if a cell is "viable" or not
#
# @param c Cell to determine if viable
#/
static inline int cellIsViable(struct Cell *const c)
{
   return ( c->generation > 2 ) ? 1 : 0; //&& c->children ) ? 1 : 0;
}

# Code reuse is good! #
static inline void incrementAndWrapPtrs( uintptr_t *wordPtr, uintptr_t *shiftPtr ) {
   if ((*shiftPtr += 4) >= SYSWORD_BITS) {
      if (++*wordPtr >= POND_DEPTH_SYSWORDS) {
	 *wordPtr = EXEC_START_WORD;
	 *shiftPtr = EXEC_START_BIT;
      } else *shiftPtr = 0;
   }
}

#*
# Output a line of comma-seperated statistics data
#
# @param clock Current clock
#/
static void doReport(FILE *file, const uint64_t clockTicks)
{
  static uint64_t lastTotalViableReplicators = 0;
  
  uintptr_t x,y;
  
  uint64_t totalActiveCells = 0;
  uint64_t totalEnergy = 0;
  uint64_t totalViableReplicators = 0;
  uintptr_t maxGeneration = 0;
  uintptr_t maxChildren = 0;
  
  for(x=0;x<POND_SIZE_X;++x) {
    for(y=0;y<POND_SIZE_Y;++y) {
      struct Cell *const c = &pond[x][y];
      if (c->energy) {
        ++totalActiveCells;
        totalEnergy += (uint64_t)c->energy;
        if ( cellIsViable( c ) )
          ++totalViableReplicators;
        if (c->generation > maxGeneration)
          maxGeneration = c->generation;
        if (c->children > maxChildren)
          maxChildren = c->children;
      }
    }
  }
  
  # Look here to get the columns in the CSV output #
  
  # The first five are here and are self-explanatory #
  fprintf(file, "%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu",
    (uint64_t)clockTicks,
    (uint64_t)totalEnergy,
    (uint64_t)totalActiveCells,
    (uint64_t)totalViableReplicators,
    (uint64_t)maxGeneration,
    (uint64_t)maxChildren,
    (uint64_t)statCounters.viableCellsReplaced,
    (uint64_t)statCounters.viableCellsKilled,
    (uint64_t)statCounters.viableCellShares,
    (uint64_t)statCounters.bababooms
    );
  
  # The next 16 are the average frequencies of execution for each
   * instruction per cell execution. #
  double totalMetabolism = 0.0;
  for(x=0;x<16;++x) {
    totalMetabolism += statCounters.instructionExecutions[x];
    fprintf(file,",%.4f",(statCounters.cellExecutions > 0.0) ? (statCounters.instructionExecutions[x] / statCounters.cellExecutions) : 0.0);
  }
  
  # The last column is the average metabolism per cell execution #
  fprintf(file,",%.4f\n",(statCounters.cellExecutions > 0.0) ? (totalMetabolism / statCounters.cellExecutions) : 0.0);
  fflush(file#stdout#);
  
  if ((lastTotalViableReplicators > 0)&&(totalViableReplicators == 0))
    fprintf(stderr,"[EVENT] Viable replicators have gone extinct. Please reserve a moment of silence.\n");
  else if ((lastTotalViableReplicators == 0)&&(totalViableReplicators > 0))
    fprintf(stderr,"[EVENT] Viable replicators have appeared!\n");
  
  lastTotalViableReplicators = totalViableReplicators;
  
  # Reset per-report stat counters #
  bzero( &statCounters, sizeof( statCounters ) );
}

#ifdef NUM_BARRIERS
# Test if cell at x,y is on a barrier, if so return 1 if on a vert. barrier, 2 if
   on a horiz. barrier. #
static inline int testBarrier( uintptr_t x, uintptr_t y ) {
   uintptr_t i;
   if ( n_barriers == 0 ) return( 0 );
   for (i=0;i<n_barriers;++i) if ( x == barriers_x[ i ] ) return( 1 );
   for (i=0;i<n_barriers;++i) if ( y == barriers_y[ i ] ) return( 2 );
   return( 0 );
}
#endif

#*
# Gets the color that a cell should be
#
# @param c Cell to get color for
# @return 8-bit color value
#/
#ifdef USE_SDL 
static inline uint8_t getColor(struct Cell *c, uintptr_t x, uintptr_t y)
{
  uintptr_t i,j,word,sum,opcode;

#ifdef NUM_BARRIERS
  if ( testBarrier( x, y ) ) return( 0xf );
#endif
  
  if (! colorOnlyViable || c->energy) {
    switch(colorScheme) {
      case KINSHIP:
        #
         * Kinship color scheme by Christoph Groth
         *
         * For cells of generation > 1, saturation and value are set to maximum.
         * Hue is a hash-value with the property that related genomes will have
         * similar hue (but of course, as this is a hash function, totally
         * different genomes can also have a similar or even the same hue).
         * Therefore the difference in hue should to some extent reflect the grade
         * of "kinship" of two cells.
         #
       if ( cellIsViable( c ) ) {
          sum = 0;
          for(i=0;i<4#POND_DEPTH_SYSWORDS#&&(c->genome[i] != ~((uintptr_t)0));++i) {
            word = c->genome[i];
            for(j=0;j<SYSWORD_BITS/4;++j,word >>= 4) {
              # We ignore 0xf's here, because otherwise very similar genomes
               * might get quite different hash values in the case when one of
               * the genomes is slightly longer and uses one more maschine
               * word. #
	       if ( i == 0 && j == 0 ) continue;
	       opcode = word & 0xf;
	       if (opcode != 0xf)
		  sum += opcode;
            }
          }
          # For the hash-value use a wrapped around sum of the sum of all
           * commands and the length of the genome. #
          return (uint8_t)((sum % 192) + 64);
        }
        return 0;
      case LOGO:
        #
         * Viable cells are color-coded by their logo.
         #
	return ( cellIsViable( c ) ) ?
	   ((uint8_t)(c->genome[0] | c->genome[0]<<4) | (uint8_t)1) : 0;
      case LINEAGE:
        #
         * Viable cells are color-coded by lineage.
         #
	return ( cellIsViable( c ) ) ?
	   (((uint8_t)c->lineage) | (uint8_t)1) : 0;
#ifdef USE_EXTRAS
      case ENERGY:
        #
         * Plot the energy available in each cell as a grey-scale.
         #
	return (((uint8_t)c->energy) | (uint8_t)1);
#endif
      case MAX_COLOR_SCHEME:
        # ... never used... to make compiler shut up. #
        break;
    }
  }
  return 0; # Cells with no energy are black #
}
#endif # USE_SDL #

static struct Cell *findClosestCell( uintptr_t *xx, uintptr_t *yy ) {
   struct Cell *cell = &pond[*xx][*yy];
   int i;
   uintptr_t x = *xx, y = *yy;
   # Find cell closest to this one if this one is not alive (for mouse presses, mostly) #
   if ( ! ( cellIsViable( cell ) ) ) {
      i = 0;
      while( i ++ < 5 && ! ( cellIsViable( cell ) ) ) {
	 cell = &pond[x+i][y]; if ( cellIsViable( cell ) ) { *xx=x+i; *yy=y; return cell; }
	 cell = &pond[x-i][y]; if ( cellIsViable( cell ) ) { *xx=x-1; *yy=y; return cell; }
	 cell = &pond[x][y+i]; if ( cellIsViable( cell ) ) { *xx=x; *yy=y+i; return cell; }
	 cell = &pond[x][y-i]; if ( cellIsViable( cell ) ) { *xx=x; *yy=y-i; return cell; }
	 cell = &pond[x+i][y+i]; if ( cellIsViable( cell ) ) { *xx=x+i; *yy=y+i; return cell; }
	 cell = &pond[x-i][y-i]; if ( cellIsViable( cell ) ) { *xx=x-i; *yy=y-i; return cell; }
	 cell = &pond[x-i][y+i]; if ( cellIsViable( cell ) ) { *xx=x-i; *yy=y+i; return cell; }
	 cell = &pond[x+i][y-i]; if ( cellIsViable( cell ) ) { *xx=x+i; *yy=y-i; return cell; }
      }
   }
   return cell;
}   

static int fullDumpCell(FILE *file, struct Cell *cell, uintptr_t x, uintptr_t y, 
			int fullDump, struct VMStruct *VM )
{
  uintptr_t wordPtr,shiftPtr,inst,stopCount,i;
  uint8_t color;

#ifdef USE_EXTRAS
  char bufferStr[POND_DEPTH_SYSWORDS], tmpStr[ 5 ]; bufferStr[ 0 ] = '\0';
#ifdef NUM_BARRIERS
  if ( testBarrier( x, y ) ) return( 0 );
#endif
#endif  
    
  if ( cellIsViable( cell ) ) {
#ifdef USE_SDL 
    color = getColor(cell,x,y);
#else
    color = cell->genome[0];
#endif
#ifdef USE_EXTRAS
    if ( fullDump ) {
       fprintf(file,
	       "ID: %llu\tparentID: %llu\tLineage: %llu\tGen: %llu\nChild: %llu\tLogo: 0x%1x\nEn: %llu\tIntEn: %llu\n",
	       (uint64_t)cell->ID, (uint64_t)cell->parentID, (uint64_t)cell->lineage,
	       (uint64_t)cell->generation, (uint64_t)cell->children, cell->genome[0]&0xf,
	       (uint64_t)cell->energy);//, (uint64_t)cell->internalEnergy);
    } else {
#endif
       fprintf(file,"%d,%d,%llu,%llu,%llu,%llu,%llu,%llu,%x,",
	       x, y, (uint64_t)cell->ID, (uint64_t)cell->parentID, (uint64_t)cell->lineage,
	       (uint64_t)cell->generation, (uint64_t)cell->children, (uint64_t)cell->energy, color);
	       
#ifdef USE_EXTRAS
    }
#endif

    wordPtr = 0;
    shiftPtr = 0;
    stopCount = 0;
    for(i=0;i<POND_DEPTH;++i) {
      inst = (cell->genome[wordPtr] >> shiftPtr) & 0xf;
      # Four STOP instructions in a row is considered the end.
       * The probability of this being wrong is *very* small, and
       * could only occur if you had four STOPs in a row inside
       * a LOOP/REP pair that's always false. In any case, this
       * would always result in our *underestimating* the size of
       * the genome and would never result in an overestimation. #
#ifdef USE_EXTRAS 
      if ( fullDump ) {
	 if ( ! ( wordPtr == 0 && shiftPtr == 0 ) )
	    fprintf(file, "\n%5d %5d 0x%x   %8s", wordPtr, shiftPtr, inst, instructionNames[ inst ] );
	 sprintf( tmpStr, "%x", (VM->outputBuf[wordPtr] >> shiftPtr) & 0xf );
	 strcat( bufferStr, tmpStr );
			     
      } else {
#endif
	 fprintf(file,"%x",inst);
#ifdef USE_EXTRAS 
      }
#endif
      if (inst == 0xf) { # STOP #
        if (++stopCount >= 4)
          break;
      } else stopCount = 0;
      incrementAndWrapPtrs( &wordPtr, &shiftPtr );
    }
    fprintf(file,"\n");
#ifdef USE_EXTRAS 
    if ( fullDump ) {
       if ( ! ( bufferStr[ 0 ] == 'f' && bufferStr[ 1 ] == 'f' ) )
	  fprintf(file, "OUTBUF: %s\n", bufferStr );
       if ( VM->loopStackPtr > 0 ) {
	  fprintf(file, "LOOPSTACK WORDPTR:  " );
	  for ( i=0; i<VM->loopStackPtr; i++ )
	     fprintf(file, "%3d ", VM->loopStack_wordPtr[i]);
	  fprintf(file, "\nLOOPSTACK SHIFTPTR: " );
	  for ( i=0; i<VM->loopStackPtr; i++ )
	     fprintf(file, "%3d ", VM->loopStack_shiftPtr[i]);
          fprintf(file, "\n");
       }
    }
#endif
    return( 1 );
  }
  return( 0 );
}

#*
# Dumps the genome of a cell to a file.
#
# @param file Destination
# @param cell Source
# @param x X coord of cell
# @param y Y coord of cell
#/
static inline int dumpCell(FILE *file, struct Cell *cell, uintptr_t x, uintptr_t y)
{
   return fullDumpCell( file, cell, x, y, 0, (struct VMStruct *) NULL );
}

#*
# Dumps all viable (generation > 2) cells to a file called <clock>.dump
#
# @param clock Clock value
#/
static int doDump(const uint64_t clockTicks)
{
  char buf[POND_DEPTH*2];
  FILE *d;
  uintptr_t x,y,wordPtr,shiftPtr,inst,i;
  int ok;
  struct Cell *pptr;

  if ( clockTicks == 0 ) return;
  
  sprintf(buf,"output/%012llu.dump.csv",clockTicks);
  d = fopen(buf,"w");
  if (!d) {
    fprintf(stderr,"[WARNING] Could not open %s for writing.\n",buf);
    return;
  }
  
  ok = 0;
  for(x=0;x<POND_SIZE_X;++x)
    for(y=0;y<POND_SIZE_Y;++y) {
#ifdef NUM_BARRIERS
       if ( testBarrier( x, y ) ) continue;
#endif
       if (dumpCell(d, &pond[x][y], x, y)) ok = 1;
    }
  
  fclose(d);
  if ( !ok ) remove( buf );
  else fprintf(stderr,"[INFO] Dumping viable cells to %s\n",buf);
  return( ok );
}

#*
# Save a screen shot to a bitmap (bmp) file.
#
# @param clock Clock value
# @param screen SDL screen pointer
#/
#ifdef USE_EXTRAS
#ifdef USE_SDL

# Store the SDL screen and pixel pitch var. globally #
static SDL_Surface *screen;
static uintptr_t sdlPitch;

static void doScreenshot(const uint64_t clockTicks)
{
  char buf[POND_DEPTH*2];

  if ( clockTicks == 0 ) return;
  
  sprintf(buf,"output/%012llu.scrn.bmp",clockTicks);
  fprintf(stderr,"[INFO] Saving screenshot to %s\n",buf);

  SDL_SaveBMP(screen, buf);
}

static void doSaveOpen(const int8_t saveOpen, uint64_t *clockTicks, uint64_t *cellIdCounter, const char *fname )
{
  char buf[POND_DEPTH*2];
  int i;
  FILE *d;
  
  if ( saveOpen == 's' ) {
     if ( *clockTicks == 0 ) return;
     sprintf(buf,"output/%012llu.dump.bin",*clockTicks);
     d = fopen(buf,"w");
     if (!d) {
       fprintf(stderr,"[WARNING] Could not open %s for writing.\n",buf);
       return;
     }
     fprintf(stderr,"[INFO] Saving memory dump to %s\n",buf);
     fwrite( clockTicks, sizeof(uint64_t), 1, d );
     fwrite( cellIdCounter, sizeof(uint64_t), 1, d );
     fwrite( &statCounters, sizeof(struct PerReportStatCounters), 1, d );
#ifdef NUM_BARRIERS
     fwrite( &n_barriers, sizeof(uintptr_t), 1, d );
     fwrite( &barriers_x, sizeof(uintptr_t), NUM_BARRIERS, d );
     fwrite( &barriers_y, sizeof(uintptr_t), NUM_BARRIERS, d );
#endif
     for (i=0; i<POND_SIZE_X; i++ )
	fwrite( &(pond[i][0]), sizeof(struct Cell), POND_SIZE_Y, d );
  } else if ( saveOpen == 'o' ) {
     d = fopen(fname,"r");
     if (!d) {
       fprintf(stderr,"[WARNING] Could not open %s for reading.\n",fname);
       return;
     }
     fprintf(stderr,"[INFO] Reading memory dump from %s\n",fname);
     fread( clockTicks, sizeof(uint64_t), 1, d );
     fread( cellIdCounter, sizeof(uint64_t), 1, d );
     fread( &statCounters, sizeof(struct PerReportStatCounters), 1, d );
#ifdef NUM_BARRIERS
     fread( &n_barriers, sizeof(uintptr_t), 1, d );
     fread( &barriers_x, sizeof(uintptr_t), NUM_BARRIERS, d );
     fread( &barriers_y, sizeof(uintptr_t), NUM_BARRIERS, d );
#endif
     for (i=0; i<POND_SIZE_X; i++ )
	fread( &(pond[i][0]), sizeof(struct Cell), POND_SIZE_Y, d );
  }
  
  fclose(d);
}
#endif
#endif

#
# Move resetting into its own function so it can be called (e.g. by 'r' key press)
#/
static void reset()
{
   uintptr_t i,x,y;

   for(x=0;x<POND_SIZE_X;++x) {
      for(y=0;y<POND_SIZE_Y;++y) {
	 pond[x][y].ID = pond[x][y].parentID = pond[x][y].lineage = 0;
	 pond[x][y].generation = pond[x][y].children = pond[x][y].energy = 0;
	 memset( pond[x][y].genome, POND_DEPTH_BYTES, 0xff );
      }
   }
  # Reset per-report stat counters #
   bzero( &statCounters, sizeof( statCounters ) );

#ifdef NUM_BARRIERS
   bzero( &barriers_x, sizeof( barriers_x ) );
   bzero( &barriers_y, sizeof( barriers_y ) );
#endif
}   

#*
# Get a neighbor in the pond
#
# @param x Starting X position
# @param y Starting Y position
# @param dir Direction to get neighbor from
# @return Pointer to neighboring cell
#/
static inline struct Cell *getNeighbor(const uintptr_t x,const uintptr_t y,const uintptr_t dir)
{
   uintptr_t xnew = x, ynew = y;
  # Space is toroidal; it wraps at edges #
  switch(dir) {
  case N_LEFT:
     xnew = (x) ? x-1 : POND_SIZE_X-1;
     break;
  case N_RIGHT:
     xnew = (x < (POND_SIZE_X-1)) ? x+1 : 0;
     break;
  case N_UP:
     ynew = (y) ? y-1 : POND_SIZE_Y-1;
     break;
  case N_DOWN:
     ynew = (y < (POND_SIZE_Y-1)) ? y+1 : 0;
     break;
  }
#ifdef NUM_BARRIERS
  # Barriers reflect (randomly) the direction, like a shoreline or mountain range #
  if ( testBarrier( xnew, ynew ) ) return( getNeighbor( x, y, getRandom() % 4 ) );
#endif
  return &pond[xnew][ynew];
}

#*
# Determines if c1 is allowed to access c2
#
# @param c2 Cell being accessed
# @param c1guess c1's "guess"
# @param sense The "sense" of this interaction
# @return True or false (1 or 0)
#/
static inline int accessAllowed(struct Cell *const c2,const uintptr_t c1guess,int sense)
{
  # Access permission is more probable if they are more similar in sense 0,
   * and more probable if they are different in sense 1. Sense 0 is used for
   * "negative" interactions and sense 1 for "positive" ones. #
  return (!c2->parentID) || ( sense ?
    ((getRandom() & 0xf) >= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) :	
    ((getRandom() & 0xf) <= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) );
}

#ifdef USE_EXTRAS
static void traceCell( struct Cell *cell, uintptr_t wordPtr, uintptr_t shiftPtr,
		       uintptr_t ptr_wordPtr, uintptr_t ptr_shiftPtr, uintptr_t inst,
		       uintptr_t reg, struct VMStruct *VM ) {
   fprintf( stderr, "\n================\nPtrs: %d %d\tPPtrs: %d %d\tReg: 0x%x\tInst: 0x%x %s\n",
	    wordPtr, shiftPtr, ptr_wordPtr, ptr_shiftPtr, reg, inst, instructionNames[ inst ] );
   fullDumpCell( stderr, cell, 0, 0, 1, VM );
}
#endif

# Execute the current cell (stored in the global VM) #
static void executeCell( struct Cell *pptr, const uintptr_t x, const uintptr_t y,
			 struct VMStruct *VM, uint64_t *cellIdCounter )
{

  # Miscellaneous variables used in the loop #
  uintptr_t i,tmp,currentWord,wordPtr,shiftPtr,inst,tmp_stopCount;

  # If this is nonzero, cell execution stops. This allows us
   * to avoid the ugly use of a goto to exit the loop. :) #
  int stop;
#ifdef USE_EXTRAS
  uintptr_t stopCount;
#endif
  
  # Virtual machine memory pointer register (which
   * exists in two parts... read the code below...) #
  uintptr_t ptr_wordPtr;
  uintptr_t ptr_shiftPtr;
  # Virtual machine loop/rep stack pointer #
  uintptr_t loopStackPtr;
  # Temporary memory pointer registers, for EJECT command #
  uintptr_t tmp_wordPtr;
  uintptr_t tmp_shiftPtr;
  
  # The main "register" #
  uintptr_t reg;
  
  # If this is nonzero, we're skipping to matching REP #
  # It is incremented to track the depth of a nested set
   * of LOOP/REP pairs in false state. #
  uintptr_t falseLoopDepth;
       
  struct Cell *tmpptr;

  # Reset the state of the VM prior to execution #
  memset( VM->outputBuf, 0xff, POND_DEPTH_BYTES );
  loopStackPtr = ptr_wordPtr = ptr_shiftPtr = reg = falseLoopDepth = 0;
  wordPtr = EXEC_START_WORD;
  shiftPtr = EXEC_START_BIT;
  stop = 0; 
#ifdef USE_EXTRAS
  stopCount = 0;
#endif
  
  # We use a currentWord buffer to hold the word we're
   * currently working on.  This speeds things up a bit
   * since it eliminates a pointer dereference in the
   * inner loop. We have to be careful to refresh this
   * whenever it might have changed... take a look at
   * the code. :) #
  currentWord = pptr->genome[0];
    
  # Keep track of how many cells have been executed #
  statCounters.cellExecutions += 1.0;

  # Core execution loop #
  while (pptr->energy&&(!stop)) {
    # Get the next instruction #
    inst = (currentWord >> shiftPtr) & 0xf;
    
    # Randomly frob either the instruction or the register with a
     * probability defined by MUTATION_RATE. This introduces variation,
     * and since the variation is introduced into the state of the VM
     * it can have all manner of different effects on the end result of
     * replication: insertions, deletions, duplications of entire
     * ranges of the genome, etc. #
    if ((getRandom() & 0xffffffff) < MUTATION_RATE) {
       tmp = getRandom(); # Call getRandom() only once for speed #
       if (tmp & 0x80) # Check for the 8th bit to get random boolean #
	  inst = tmp & 0xf; # Only the first four bits are used here #
       else reg = tmp & 0xf;
    }
      
    # Each instruction processed costs one unit of energy #
    #if ( pptr->internalEnergy )
       --pptr->internalEnergy;
       else# --pptr->energy;
      
    # Execute the instruction #
    if (falseLoopDepth) {
       # Skip forward to matching REP if we're in a false loop. #
       if (inst == 0x9) # Increment false LOOP depth #
	  ++falseLoopDepth;
       else if (inst == 0xa) # Decrement on REP #
	  --falseLoopDepth;
    } else {
       # If we're not in a false LOOP/REP, execute normally #
        
       # Keep track of execution frequencies for each instruction #
       statCounters.instructionExecutions[inst] += 1.0;
        
       switch(inst) {
       case 0x0: # ZERO: Zero VM state registers #
	  reg = ptr_wordPtr = ptr_shiftPtr = 0;
	  break;
       case 0x1: # FWD: Increment the pointer (wrap at end; don't skip the LOGO -- why not?) #
	  if ((ptr_shiftPtr += 4) >= SYSWORD_BITS) {
	     if (++ptr_wordPtr >= POND_DEPTH_SYSWORDS)
		ptr_wordPtr = 0;
	     ptr_shiftPtr = 0;
	  }
	  break;
       case 0x2: # BACK: Decrement the pointer (wrap at beginning) #
	  if (ptr_shiftPtr)
	     ptr_shiftPtr -= 4;
	  else {
	     if (ptr_wordPtr)
		--ptr_wordPtr;
	     else ptr_wordPtr = POND_DEPTH_SYSWORDS - 1;
	     ptr_shiftPtr = SYSWORD_BITS - 4;
	  }
	  break;
       case 0x3: # INC: Increment the register #
	  reg = (reg + 1) & 0xf;
	  break;
       case 0x4: # DEC: Decrement the register #
	  reg = (reg - 1) & 0xf;
	  break;
       case 0x5: # READG: Read into the register from genome #
	  reg = (pptr->genome[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
	  break;
       case 0x6: # WRITEG: Write out from the register to genome #
	  pptr->genome[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
	  pptr->genome[ptr_wordPtr] |= reg << ptr_shiftPtr;
	  currentWord = pptr->genome[wordPtr]; # Must refresh in case this changed! #
	  break;
       case 0x7: # READB: Read into the register from buffer #
	  reg = (VM->outputBuf[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
	  break;
       case 0x8: # WRITEB: Write out from the register to buffer #
	  VM->outputBuf[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
	  VM->outputBuf[ptr_wordPtr] |= reg << ptr_shiftPtr;
	  break;
       case 0x9: # LOOP: Jump forward to matching REP if register is zero #
	  if (reg) {
	     if (loopStackPtr >= POND_DEPTH)
		stop = 1; # Stack overflow ends execution #
	     else {
		VM->loopStack_wordPtr[loopStackPtr] = wordPtr;
		VM->loopStack_shiftPtr[loopStackPtr] = shiftPtr;
		++loopStackPtr;
	     }
	  } else falseLoopDepth = 1;
	  break;
       case 0xa: # REP: Jump back to matching LOOP if register is nonzero #
	  if (loopStackPtr) {
	     --loopStackPtr;
	     if (reg) {
		wordPtr = VM->loopStack_wordPtr[loopStackPtr];
		shiftPtr = VM->loopStack_shiftPtr[loopStackPtr];
		currentWord = pptr->genome[wordPtr];
		# This ensures that the LOOP is rerun #
		continue;
	     }
	  }
	  break;
       case 0xb: # XCHG: Skip next instruction and exchange value of register with it #
	  incrementAndWrapPtrs( &wordPtr, &shiftPtr );
	  tmp = reg;
	  reg = (pptr->genome[wordPtr] >> shiftPtr) & 0xf;
	  pptr->genome[wordPtr] &= ~(((uintptr_t)0xf) << shiftPtr);
	  pptr->genome[wordPtr] |= tmp << shiftPtr;
	  currentWord = pptr->genome[wordPtr];
	  break;
       case 0xc: # EJECT: eject output buffer into neighbor at location 0.
	   * Get neighbor direction from register and access permission from first
	   * value in buffer. As with old KILL, there is a penalty on failure. #
	  if ((VM->outputBuf[0] & 0xff) != 0xff) {
	     tmpptr = getNeighbor(x,y,reg&3);
	     if (accessAllowed(tmpptr,VM->outputBuf[0]&0xf,0)) {
		# Log it if we're replacing a viable cell #
		if ( cellIsViable( tmpptr ) )
		   ++statCounters.viableCellsReplaced;
        
		tmpptr->ID = ++ *cellIdCounter;
		tmpptr->parentID = pptr->ID;
		tmpptr->generation = pptr->generation + 1;
		tmpptr->children = 0;
		pptr->children = pptr->children + 1;
#ifndef USE_EXTRAS
		memcpy( tcmpptr->genome, VM->outputBuf, POND_DEPTH_BYTES );
		tmpptr->lineage = pptr->lineage; # Lineage is copied in offspring #
#else
		# NOTE both these below are WRONG - we want to copy starting from the 2nd instr in
		   the genome, even for the first case -- because the first instr. is used for the
		   location of injection! #
		if ( ptr_wordPtr == 0 && ptr_shiftPtr == 0 ) {
		   # Insertion location is at zero - it's a KILL #
		   memcpy( tmpptr->genome, VM->outputBuf, POND_DEPTH_BYTES );
		   tmpptr->lineage = pptr->lineage; # Lineage only follows the logo (instr. 1) #
		} else {
		   # Insertion location is at non-zero - it's sex! #
		   # Copy output buffer into neighbor, at location given by ptr_wordPtr, ptr_shiftPtr #
		   if ( cellIsViable( tmpptr ) ) # Log it if we're having sex #
		      ++statCounters.bababooms;
		   tmp_stopCount = 0;
		   tmp_wordPtr = 0; //EXEC_START_WORD;
		   tmp_shiftPtr = 0; //EXEC_START_SHIFT;
		   
		   for(i=0;i<POND_DEPTH;++i) {
		      inst = (VM->outputBuf[tmp_wordPtr] >> tmp_shiftPtr) & 0xf;
		   
		      tmpptr->genome[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
		      tmpptr->genome[ptr_wordPtr] |= inst << ptr_shiftPtr;
	
		      if (inst == 0xf) { // STOP 
			 if (++tmp_stopCount >= 4) 
			    break;
		      } else tmp_stopCount = 0;

		      # Treat genome as circular - when copying stops at end, it continues at beginning. #
		      incrementAndWrapPtrs( &ptr_wordPtr, &ptr_shiftPtr );
		      incrementAndWrapPtrs( &tmp_wordPtr, &tmp_shiftPtr );
		   }
		}
#endif
		// Reset the VM state 
		reg = ptr_wordPtr = ptr_shiftPtr = 0;
		memset( VM->outputBuf, 0xff, POND_DEPTH_BYTES );
	     }
#ifdef USE_EXTRAS
	     else if ( cellIsViable( tmpptr ) ) {
		tmp = pptr->energy / FAILED_EJECT_PENALTY;
		if (pptr->energy > tmp)
		   pptr->energy -= tmp;
		else pptr->energy = 0;
	     }
#endif
	  }
	  break;
       case 0xd: # Metabolize energy? (Make energy available to use as internalEnergy) #
#ifdef USE_EXTRAS
	  #if ( pptr->energy >= METABOLITE_AMOUNT ) {
	     pptr->internalEnergy += METABOLITE_AMOUNT;
	     pptr->energy -= METABOLITE_AMOUNT;
	  } else if ( pptr->energy > 0 ) {
	     pptr->internalEnergy += pptr->energy;
	     pptr->energy = 0;
	  }
	  #
	  break; # KILL Not needed anymore because of EJECT? #
#endif
       case 0xe: # SHARE: Equalize energy between self and neighbor if allowed.
		  Note this also can be "stealing". #
	  tmpptr = getNeighbor(x,y,reg&3);
	  if ( ! tmpptr->children || tmpptr->parentID ||
	       (pptr->energy > pptr->energy && accessAllowed(tmpptr,reg,1)) || # Sharing #
	       (pptr->energy < pptr->energy && accessAllowed(tmpptr,reg,0))) { # Stealing #

	     if ( cellIsViable( tmpptr ) )
		++statCounters.viableCellShares;

	     tmp = pptr->energy + tmpptr->energy;
	     tmpptr->energy = tmp / 2;
	     pptr->energy = tmp - tmpptr->energy;
	  }
	  break;
       case 0xf: # STOP: End execution #
#ifdef USE_EXTRAS
	  if (++stopCount >= 4) 
#endif
	  stop = 1;
	  break;
       }
#ifdef USE_EXTRAS
       if ( inst != 0xf )
	  stopCount = 0;
#endif
    }

#ifdef USE_EXTRAS
    if ( pptr == tracePtr && tracePtr != (struct Cell *) NULL && cellIsViable( tracePtr ) ) {
       VM->loopStackPtr = loopStackPtr;
       traceCell( tracePtr, wordPtr, shiftPtr, ptr_wordPtr, ptr_shiftPtr, inst, reg, VM );
    }
#endif
    
    # Advance the shift and word pointers, and loop around
     * to the beginning at the end of the genome. #
    incrementAndWrapPtrs( &wordPtr, &shiftPtr );
    currentWord = pptr->genome[wordPtr];
  }
  # Do we only want to print out one run-through of the traceCell's program? #
  if ( pptr == tracePtr && tracePtr != (struct Cell *) NULL ) tracePtr = (struct Cell *)NULL;
}

# Introduce a random cell or mutation somewhere with a given energy level #
# This is called seeding, and introduces both energy and
# entropy into the substrate. This happens every INFLOW_FREQUENCY
# clock ticks. #
static void mutateCell( struct Cell *pptr, const uintptr_t x, const uintptr_t y,
			 struct VMStruct *VM, uint64_t *cellIdCounter ) {
   uintptr_t i;
#ifdef USE_EXTRAS   
   if ( cellIsViable(pptr) && (getRandom() & 0xffffffff) < MUTATION_RATE ) {
      # Random mutation in the genome of this living cell (like UV radiation) #
      pptr->genome[getRandom() % POND_DEPTH_SYSWORDS] = getRandom();
   } else {
#endif
      for(i=0;i<POND_DEPTH_SYSWORDS;++i) 
	 pptr->genome[i] = getRandom();
      pptr->ID = *cellIdCounter;
      pptr->parentID = 0;
      pptr->lineage = *cellIdCounter;
      pptr->generation = pptr->children = 0;
      ++ *cellIdCounter;
#ifdef USE_EXTRAS   
   }
#endif
}	

#ifdef USE_SDL
static int processSDLEvents() {
  uintptr_t i,x,y;
  int tmp;
  struct Cell *tmpptr;
  SDL_Event sdlEvent;

  while (SDL_PollEvent(&sdlEvent)) {
     if (sdlEvent.type == SDL_QUIT) {
	fprintf(stderr,"[QUIT] Quit signal received!\n");
	return(1);
     } else if (sdlEvent.type == SDL_MOUSEBUTTONDOWN) {
	switch (sdlEvent.button.button) {
	case SDL_BUTTON_LEFT:
	   x = sdlEvent.button.x; y = sdlEvent.button.y;
	   fprintf(stderr,"[INTERFACE] Genome of cell at (%d, %d):\n",x, y );
	   if ( x < POND_SIZE_X && x > 0 && y < POND_SIZE_Y && y > 0 ) {
	      tmpptr = findClosestCell( &x, &y );
	      if ( tmpptr != (struct Cell *) NULL )
		 dumpCell(stderr, tmpptr, x, y);
	   }
	   break;
	}
     } else if ( sdlEvent.type == SDL_KEYDOWN ) {
	switch( sdlEvent.key.keysym.sym ) {
	case 'q': 
	   fprintf(stderr,"[QUIT] Quit signal received!\n");
	   return(1);
	   break;
	case 'p':
	   fprintf(stderr,"[PAUSE] Pause signal received!\n");
	   paused = !paused;
	   break;
	case '1':
	case '2':
	case '3':
	case '4':
	   colorScheme = sdlEvent.key.keysym.sym - '1';
	   fprintf(stderr,"[COLOR SCHEME] Switching to color scheme \"%s\".\n",
		   colorSchemeName[colorScheme]);
	   for (y=0;y<POND_SIZE_Y;++y) {
	      for (x=0;x<POND_SIZE_X;++x)
		 ((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = getColor(&pond[x][y],x,y);
	   }
	   break;
#ifdef USE_EXTRAS
	case 'r':
	   fprintf(stderr,"[RESTART] Restart signal received!\n");
	   reset();
	   break;
	case 'd':
	   fprintf(stderr,"[DUMP] Dump signal received!\n");
	   if (doDump(clockTicks))
	      doScreenshot(clockTicks);
	   doReport(stderr,clockTicks);
	   break;
	case 's':
	   fprintf(stderr,"[SAVE] Save signal received!\n");
	   doSaveOpen( 's', &clockTicks, &cellIdCounter, (char *)NULL );
	   break;
	case '-':
	   if ( adjustedMeanEnergy > -INFLOW_RATE_BASE + 101 )
	      adjustedMeanEnergy -= 100;
	   else adjustedMeanEnergy -= 10;
	   if (adjustedMeanEnergy <= -INFLOW_RATE_BASE)
	      adjustedMeanEnergy = -INFLOW_RATE_BASE+1;
	   fprintf(stderr,"[ENERGY] Mean energy input decreased to %d!\n",
		   INFLOW_RATE_BASE+adjustedMeanEnergy);
	   break;
	case '=':
	   if ( adjustedMeanEnergy >= -INFLOW_RATE_BASE + 100 )
	      adjustedMeanEnergy += 100;
	   else adjustedMeanEnergy += 10;
	   fprintf(stderr,"[ENERGY] Mean energy input increased to %d!\n",
		   adjustedMeanEnergy+INFLOW_RATE_BASE);
	   break;
	case 'a':
	   colorOnlyViable = ! colorOnlyViable;
	   fprintf(stderr,"[COLOR SCHEME] Coloring only cells that are viable: \"%d\".\n",
		   colorOnlyViable);
	   break;
	case 'i': # Insert a genome of your choosing at the mouse location #
	   SDL_GetMouseState((int*) &x, (int*) &y);
	   if ( x >= POND_SIZE_X || x <= 0 || y >= POND_SIZE_Y && y <= 0 ) break;
	   fprintf( stderr,"[INSERT] Enter the text of a genome to insert here (%d, %d):\n", x, y );
	   struct Cell *pptr = &pond[x][y];
	   for(i=0;i<POND_DEPTH_SYSWORDS;++i) 
	      pptr->genome[i] = getRandom();
	   pptr->ID = cellIdCounter;
	   pptr->parentID = 0;
	   pptr->lineage = cellIdCounter;
	   pptr->generation = pptr->children = 10;
	   ++ cellIdCounter;
	   char tmp_c[ 5 ];
	   uintptr_t inst, wordPtr = EXEC_START_WORD, shiftPtr = EXEC_START_BIT;
	   i = 1;
	   while( ( tmp = getchar() ) != '\n' ) {
	      sprintf( tmp_c, "0x%c", tmp );
	      inst = strtoq( tmp_c, (char **) NULL, 16 ) & 0xf;
	      //fprintf( stderr,"INPUT: %d %c %s %x %d %d\n",i,i,tmp_c,inst,wordPtr,shiftPtr);
	      pptr->genome[wordPtr] &= ~(((uintptr_t)0xf) << shiftPtr);
	      pptr->genome[wordPtr] |= inst << shiftPtr;
	      incrementAndWrapPtrs( &wordPtr, &shiftPtr );
	   }	
	   fprintf(stderr,"[INTERFACE] Genome of cell at (%d, %d):\n", x, y );
	   dumpCell( stderr, pptr, (uintptr_t) x, (uintptr_t) y );
	   pptr->generation = pptr->children = 0;
	   break;
	case 't': # Trace the current bug (to stderr) until 't' is pressed again #
	   SDL_GetMouseState((int*) &x, (int*) &y);
	   if ( x >= POND_SIZE_X || x <= 0 || y >= POND_SIZE_Y && y <= 0 ) break;
	   if ( x >= POND_SIZE_X || y >= POND_SIZE_Y ) break;
	   if ( tracePtr == &pond[x][y] ) {
	      fprintf(stderr,"[INTERFACE] Stopped tracing cell.\n" );
	      tracePtr = (struct Cell *)NULL;
	   } else {
	      fprintf(stderr,"[INTERFACE] Tracing cell at (%d, %d):\n", x, y );
	      tracePtr = findClosestCell( &x, &y );
	   }
	   break;
#ifdef NUM_BARRIERS
	case 'b':
	   if ( ! SDL_GetModState() ) { # Lower-case 'b' #
	      if ( n_barriers < NUM_BARRIERS ) {
		 SDL_GetMouseState((int*) &x, (int*) &y);
		 barriers_x[ n_barriers ] = x; //getRandom() % POND_SIZE_X;
		 barriers_y[ n_barriers ++ ] = y; //getRandom() % POND_SIZE_Y;
		 fprintf(stderr,"[BARRIERS] Adding barriers at (%d, %d): %d.\n", x, y, n_barriers );
	      }
	   } else { # Capital 'B' #
	      if ( n_barriers > 0 ) {
		 barriers_x[ -- n_barriers ] = barriers_y[ n_barriers ] = 0;
		 fprintf(stderr,"[BARRIERS] Removing barriers: %d.\n", n_barriers );
	      }
	   }
	   break;
#endif
#endif
	}
     }
  }  
}

static inline void updateSDLPixels( int x, int y, int neighbors ) {
   if (SDL_MUSTLOCK(screen))
      SDL_LockSurface(screen);
   ((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = getColor(&pond[x][y],x,y);
   if ( neighbors ) {
      if (x) {
	 ((uint8_t *)screen->pixels)[(x-1) + (y * sdlPitch)] = getColor(&pond[x-1][y],x-1,y);
	 if (x < (POND_SIZE_X-1))
	    ((uint8_t *)screen->pixels)[(x+1) + (y * sdlPitch)] = getColor(&pond[x+1][y],x+1,y);
	 else ((uint8_t *)screen->pixels)[y * sdlPitch] = getColor(&pond[0][y],0,y);
      } else {
	 ((uint8_t *)screen->pixels)[(POND_SIZE_X-1) + (y * sdlPitch)] = getColor(&pond[POND_SIZE_X-1][y],POND_SIZE_X-1,y);
	 ((uint8_t *)screen->pixels)[1 + (y * sdlPitch)] = getColor(&pond[1][y],1,y);
      }
      if (y) {
	 ((uint8_t *)screen->pixels)[x + ((y-1) * sdlPitch)] = getColor(&pond[x][y-1],x,y-1);
	 if (y < (POND_SIZE_Y-1))
	    ((uint8_t *)screen->pixels)[x + ((y+1) * sdlPitch)] = getColor(&pond[x][y+1],x,y+1);
	 else ((uint8_t *)screen->pixels)[x] = getColor(&pond[x][0],x,0);
      } else {
	 ((uint8_t *)screen->pixels)[x + ((POND_SIZE_Y-1) * sdlPitch)] = getColor(&pond[x][POND_SIZE_Y-1],x,POND_SIZE_Y-1);
	 ((uint8_t *)screen->pixels)[x + sdlPitch] = getColor(&pond[x][1],x,1);
      }
   }
   if (SDL_MUSTLOCK(screen))
      SDL_UnlockSurface(screen);
}
#endif

int doMainLoop( void *data ) {
  uintptr_t i,x,y;
  int tmp;
  struct Cell *tmpptr;

  # Global VM #
  struct VMStruct VM;

  # Main loop #
  for(;;) {
#ifdef USE_SDL
    # SDL display is refreshed every SDL_UPDATE_FREQUENCY #
    if (!(clockTicks % SDL_UPDATE_FREQUENCY)) {
      if ( processSDLEvents() ) # Quit was chosen by user if returns 1 #
	 return( 0 );
      
      if (!paused)
	 SDL_UpdateRect(screen,0,0,POND_SIZE_X,POND_SIZE_Y);
    }
    if (paused) {
       sleep(2);
       continue;
    }       
#endif # USE_SDL #
    
    # Stop at STOP_AT if defined #
#ifdef STOP_AT
    if (clockTicks >= STOP_AT) {
      # Also do a final dump if dumps are enabled #
#ifdef DUMP_FREQUENCY
      doDump(clockTicks);
#ifdef USE_EXTRAS
#ifdef USE_SDL
      doScreenshot(clockTicks);
#endif
#endif
#endif # DUMP_FREQUENCY #
      fprintf(stderr,"[QUIT] STOP_AT clock value reached\n");
      break;
    }
#endif # STOP_AT #

    # Increment clock and run reports periodically #
    # Clock is incremented at the start, so it starts at 1 #
    if (!(++clockTicks % REPORT_FREQUENCY))
      doReport(stdout,clockTicks);

    # Periodically dump the viable population if defined #
#ifdef DUMP_FREQUENCY
    if (!(clockTicks % DUMP_FREQUENCY)) {
      tmp = doDump(clockTicks);
#ifdef USE_EXTRAS
#ifdef USE_SDL
      if (tmp) doScreenshot(clockTicks);
#endif
#endif
    }
#endif # DUMP_FREQUENCY #
    
    # Pick a random cell to add a new random genome to or to mutate #
    if (!(clockTicks % (INFLOW_FREQUENCY))) {
      x = getRandom() % POND_SIZE_X;
      y = getRandom() % POND_SIZE_Y;
#ifdef NUM_BARRIERS
      if ( ! testBarrier( x, y ) ) {
#endif
	 tmpptr = &pond[x][y];
	 mutateCell( tmpptr, x, y, &VM, &cellIdCounter );
	 # Add energy too #
#ifdef INFLOW_RATE_VARIATION
	 tmpptr->energy += INFLOW_RATE_BASE + (getRandom() % INFLOW_RATE_VARIATION);
#else
	 tmpptr->energy += INFLOW_RATE_BASE;
#endif # INFLOW_RATE_VARIATION #
#ifdef USE_EXTRAS
	 tmpptr->energy += adjustedMeanEnergy;
#endif
	 # Update the random cell on SDL screen if viz is enabled #
#ifdef USE_SDL
	 updateSDLPixels( x, y, 0 );
#endif # USE_SDL #
#ifdef NUM_BARRIERS
      }
#endif
    }

    # Pick a random cell to execute #
    x = getRandom() % POND_SIZE_X;
    y = getRandom() % POND_SIZE_Y;
#ifdef NUM_BARRIERS
    if ( testBarrier( x, y ) ) continue;
#endif
    tmpptr = &pond[x][y];

    # Execute the cell #
    executeCell( tmpptr, x, y, &VM, &cellIdCounter );

    # Update the neighborhood on SDL screen to show any changes. #
#ifdef USE_SDL
    updateSDLPixels( x, y, 1 );
#endif # USE_SDL #
  }
  return( 0 );
}   

#*
# Main method
#
# @param argc Number of args
# @param argv Argument array
#/
int main(int argc,char **argv)
{
  uintptr_t i;
   
  # Seed and init the random number generator #
  init_genrand(time(NULL));
  for(i=0;i<1024;++i)
    getRandom();

  # Clear the pond and initialize all genomes
   * to 0xffff and all counters to 0 ... #
  reset();
  
  # Set up SDL if we're using it #
#ifdef USE_SDL
  if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
    fprintf(stderr,"*** Unable to init SDL: %s ***\n",SDL_GetError());
    exit(1);
  }
  atexit(SDL_Quit);
  SDL_WM_SetCaption("nanopond","nanopond");
  screen = SDL_SetVideoMode(POND_SIZE_X,POND_SIZE_Y,8,SDL_SWSURFACE);
  if (!screen) {
    fprintf(stderr, "*** Unable to create SDL window: %s ***\n", SDL_GetError());
    exit(1);
  }
  sdlPitch = screen->pitch;
#endif # USE_SDL #
   
#ifdef USE_EXTRAS
#ifdef USE_SDL
  if ( argc > 1 && ! strcmp( "-f", argv[ 1 ] ) )
     doSaveOpen( 'o', &clockTicks, &cellIdCounter, argv[ 2 ] );
#endif
#endif

#ifdef N_THREADS
  int status = 0;
  SDL_Thread *thread1, *thread2;
  SDL_mutex *lock = SDL_CreateMutex();
  thread1 = SDL_CreateThread( doMainLoop, (void *) lock );
  if ( N_THREADS > 1 ) thread2 = SDL_CreateThread( doMainLoop, (void *) lock );	
  SDL_WaitThread( thread1, &status );
  if ( N_THREADS > 1 ) SDL_WaitThread( thread2, &status );
  SDL_DestroyMutex( lock );
  exit( status );
#else
  doMainLoop( (void *) "Main loop" );
#endif
  
  exit(0);
  return 0; # Make compiler shut up #
}
