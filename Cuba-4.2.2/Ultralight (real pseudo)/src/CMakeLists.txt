CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

add_executable(Rates         CalcRate.cpp          Rate.cpp Utilities.cpp CalcAmpRateRelic.cpp)

target_link_libraries(Rates ${GSL_LIBRARIES} ${Boost_LIBRARIES} cubal)

install(TARGETS Rates DESTINATION ${PROJECT_BINARY_DIR}/bin)