# ==============================================================================
# LBFGS CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(FLBFGS_OBJ
        lbfgs.f
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(flbfgs_static STATIC ${FLBFGS_OBJ})
    SET_TARGET_PROPERTIES(flbfgs_static PROPERTIES
                            OUTPUT_NAME flbfgs
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(flbfgs_shared SHARED ${FLBFGS_OBJ})
    SET_TARGET_PROPERTIES(flbfgs_shared PROPERTIES
                            OUTPUT_NAME flbfgs
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${LBFGS_VERSION})
    TARGET_LINK_LIBRARIES(flbfgs_shared
                        ${BLAS_LIB_NAME})
    INSTALL(TARGETS flbfgs_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

# ==============================================================================

# objects in library -----------------------------------------------------------
SET(FLBFGS90_OBJ
        lbfgs.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(flbfgs90_static STATIC ${FLBFGS90_OBJ})
    SET_TARGET_PROPERTIES(flbfgs90_static PROPERTIES
                            OUTPUT_NAME flbfgs90
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(flbfgs90_shared SHARED ${FLBFGS90_OBJ})
    SET_TARGET_PROPERTIES(flbfgs90_shared PROPERTIES
                            OUTPUT_NAME flbfgs90
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${LBFGS_VERSION})
    TARGET_LINK_LIBRARIES(flbfgs90_shared
                        ${BLAS_LIB_NAME})
    INSTALL(TARGETS flbfgs90_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

