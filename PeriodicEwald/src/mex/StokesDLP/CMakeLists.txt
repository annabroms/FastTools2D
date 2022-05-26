include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_stokes_dlp_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_kspace.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_pressure_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_pressure_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_pressure_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_pressure_kspace.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_gradient_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_gradient_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_gradient_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_gradient_kspace.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_pressure_grad_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_pressure_grad_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_pressure_grad_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_pressure_grad_kspace.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_vorticity_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_vorticity_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_vorticity_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_vorticity_kspace.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_stress_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_stress_real.cpp
	LINK_TO gomp
)

matlab_add_mex(
	NAME mex_stokes_dlp_stress_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_stress_kspace.cpp
	LINK_TO gomp
)

set_target_properties(mex_stokes_dlp_kspace PROPERTIES R2017b R2017b)
