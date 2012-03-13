add_executable (pross pross.cc ../base/basal_strength/basal_resistance.cc)
target_link_libraries(pross pismstressbalance pismbase)

install (TARGETS pross
  DESTINATION ${Pism_BIN_DIR})
