execs <- c("fl_statistics", "fl_generate", "fl_draw")
if(WINDOWS) execs <- paste0(execs, ".exe")
if (any(file.exists(execs))) {
  dest <- file.path(R_PACKAGE_DIR,  paste0('exec', R_ARCH))
  message("Installing ", paste(execs, collapse = " "), " to ", dest)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(execs, dest, overwrite = TRUE)
}

