# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
add_warpx_test(
    test_rz_spacecraft_charging_picmi  # name
    RZ  # dims
    2  # nprocs
    inputs_test_rz_spacecraft_charging_picmi.py  # inputs
    "analysis.py diags/diag1/"  # analysis
    "analysis_default_regression.py --path diags/diag1/"  # checksum
    OFF  # dependency
)
endif()
