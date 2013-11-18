# ============================================================================ #
#  Search-A-Gram Instagram Social Network Searcher                             #
#  http://gitlab.eugenecys.com/yjwong/cs3246_assignment3_searchagram           #
# ============================================================================ #
#  Copyright 2013 Wong Yong Jie                                                #
#                                                                              #
#  Licensed under the Apache License, Version 2.0 (the "License");             #
#  you may not use this file except in compliance with the License.            #
#  You may obtain a copy of the License at                                     #
#                                                                              #
#  http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                              #
#  Unless required by applicable law or agreed to in writing, software         #
#  distributed under the License is distributed on an "AS IS" BASIS,           #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
#  See the License for the specific language governing permissions and         #
#  limitations under the License.                                              #
# ============================================================================ #

cmake_minimum_required (VERSION 2.8.11)
project (lshkit)

# Set the cmake modules path.
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${lshkit_SOURCE_DIR}/cmake/Modules")

# Require GSL for building.
find_package (GSL)

# Add additional include directories.
include_directories (${lshkit_SOURCE_DIR}/include
  ${GSL_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR})

add_subdirectory (src)

# Require boost for demo programs.
option (BUILD_TOOLS "Build the tools for LSHKIT." OFF)
if (BUILD_TOOLS)
  find_package (Boost 1.35 COMPONENTS program_options)
  if (NOT Boost_FOUND)
    set (Boost_INCLUDE_DIR ${lshkit_SOURCE_DIR}/3rd-party/boost)
    set (Boost_LIBRARIES boost_program_options)
  endif ()
  
  add_subdirectory (tools)
endif ()

