# .show_license_message shows once and respects suppression

    Code
      .show_license_message(suppress = FALSE)
    Message
      Data from IMGT (<url> is for non-commercial use only, under a CC-BY-NC-ND 4.0 license. Attribution is required.
      HLA data from IPD-IMGT/HLA (<url> is provided under the CC-BY-ND 4.0 license.

# listIMGT() returns files or empty vector with message

    Code
      with_mocked_bindings({
        out <- listIMGT()
        expect_identical(out, character(0))
      }, .get_cache_dir = function() tmp_noexist, .package = pkg)
    Message
      Cache directory does not exist. No datasets found.

