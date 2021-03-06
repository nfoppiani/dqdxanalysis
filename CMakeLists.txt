# ======================================================================
#  dqdxanalysis main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/dqdxanalysis/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(dqdxanalysis)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  EXTRA_CXX_FLAGS -std=c++14 )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larsoft v1_00_00 )
find_ups_product( nutools v1_00_00 )
find_ups_product( art v1_09_00 )
find_ups_product( cetbuildtools v4_14_01 )


find_ups_product( larcore v1_00_00 )
find_ups_product( larcoreobj )
find_ups_product( lardata v1_00_00 )
find_ups_product( lardataobj )
find_ups_product( larevt v1_00_00 )
find_ups_product( larsim v1_00_00 )
find_ups_product( larreco v1_00_00 )
find_ups_product( larpandora v1_00_00 )
find_ups_product( larana v1_00_00 )
find_ups_product( larpandora v1_00_00 )
find_ups_product( nutools v1_00_00 )
find_ups_product( postgresql v9_1_5 )
find_ups_boost( v1_53_0 )
find_ups_product( uboonedaq_datatypes )
find_ups_product( ifdh_art )
find_ups_product( gsl )

cet_find_library( GSL NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )
cet_find_library( PQ NAMES pq PATHS ENV POSTGRESQL_LIBRARIES NO_DEFAULT_PATH )
cet_find_library( UBDut NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )


# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(dqdxanalysis)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
