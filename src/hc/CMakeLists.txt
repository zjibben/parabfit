project(Pececillo-HC)

include_directories(${Pececillo_SOURCE_DIR}/src/lib)
include_directories(${Pececillo_MODULE_DIR} ${Petaca_MODULE_DIR})

set(SOURCE_FILES
    HC_model_type.F90
    HC_precon_type.F90
    HC_AE_solver_type.F90
    HC_norm_type.F90
    HC_idaesol_model_type.F90
    HC_solver_type.F90
    HC_sim_type.F90)

# Create a library (so we can unit test)
add_library(pececillo-hc ${SOURCE_FILES})
target_link_libraries(pececillo-hc pececillo)

# HC driver
add_executable(pececillo-HC HC_main.F90)
target_link_libraries(pececillo-HC pececillo-hc)

add_subdirectory(test)

install(TARGETS pececillo-hc DESTINATION lib)
install(TARGETS pececillo-HC DESTINATION bin)
