# Install script for directory: /home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin" TYPE EXECUTABLE FILES "/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/src/Rates")
  if(EXISTS "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates"
         OLD_RPATH "/home/parracho/Desktop/Project/Cuba-4.2.2:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight_/build/bin/Rates")
    endif()
  endif()
endif()

