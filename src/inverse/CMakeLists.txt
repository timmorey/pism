# Classes for inverse problems.

set( PISM_INVERSE_SOURCES 
  ./InvTaucParameterization.cc
  ./InvSSAForwardProblem.cc
  ./InvSSAForwardProblem_dep.cc
  ./Functional.cc
  ./L2NormFunctional.cc
  ./H1NormFunctional.cc
  ./MeanSquareFunctional.cc
  ./TwoBlockVec.cc
  ./InvSSATikhonovGN.cc  
  )

if(TAO_FOUND)
  set(PISM_INVERSE_SOURCES ${PISM_INVERSE_SOURCES}
    ./TaoUtil.cc
    ./InvSSATikhonov.cc
    ./InvSSATikhonovLCL.cc
    )
endif()

add_library (pisminverse
  ${PISM_INVERSE_SOURCES}
  )

target_link_libraries (pisminverse pismbase pismutil)

install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})

# Regression tests

add_custom_target(test_inversion
  ${CMAKE_MAKE_PROGRAM} -C${CMAKE_CURRENT_BINARY_DIR} test)

macro(inversion_test name script)
  add_test("${name}:${script}" ${CMAKE_CURRENT_SOURCE_DIR}/${script}
    ${PROJECT_BINARY_DIR} ${MPIEXEC} ${PROJECT_SOURCE_DIR} -python ${PYTHON_EXECUTABLE})
endmacro()

if(Pism_BUILD_PYTHON_BINDINGS)

  execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "import siple"
    RESULT_VARIABLE IMPORT_SIPLE_ERRCODE
    OUTPUT_QUIET ERROR_QUIET)

  if (${IMPORT_SIPLE_ERRCODE} EQUAL 0) # siple is installed

    message (STATUS "Enabling regression tests that use siple")

    # copy the scripts the next test needs:
    foreach (FILE
        ${PROJECT_SOURCE_DIR}/examples/inverse/make_synth_ssa.py
        ${PROJECT_SOURCE_DIR}/examples/inverse/vel2tauc.py
        ${PROJECT_SOURCE_DIR}/examples/inverse/verify_ssa_inv.py
        test/build_tiny.py)
      get_filename_component(OUTPUT ${FILE} NAME)
      configure_file (${FILE} ${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT} COPYONLY)
    endforeach()

    inversion_test (inversion_simple test/test_tiny.sh)

  endif()
endif()
