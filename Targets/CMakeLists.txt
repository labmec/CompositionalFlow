add_executable(Compositional main_compositional.cpp
              TPZCompositional.cpp)
target_link_libraries(Compositional PRIVATE NeoPZ::pz)
