/*
 * stat_compat.c - Fortran STAT() compatibility shim for non-GNU compilers.
 *
 * gfortran provides a non-standard STAT(name, buf) extension function whose
 * runtime symbol is stat_().  LLVM Flang does not ship this symbol, so we
 * supply our own implementation that follows the same Fortran 77 calling
 * convention (trailing hidden CHARACTER length as int64_t).
 *
 * Only compiled when the Fortran compiler is not GNU (see CMakeLists.txt).
 */

#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/*
 * stat_(name, buf, name_len)
 *
 * Fortran calling convention for CHARACTER external functions:
 *   name     - pointer to Fortran CHARACTER data (not NUL-terminated)
 *   buf      - pointer to INTEGER(4) array of size >= 13
 *   name_len - hidden trailing argument: length of `name` (int64_t)
 *
 * Fills buf[] with the same 13 fields as gfortran's stat() extension:
 *   buf[0]  = st_dev,   buf[1]  = st_ino,   buf[2]  = st_mode,
 *   buf[3]  = st_nlink, buf[4]  = st_uid,   buf[5]  = st_gid,
 *   buf[6]  = st_rdev,  buf[7]  = st_size,
 *   buf[8]  = st_atime, buf[9]  = st_mtime, buf[10] = st_ctime,
 *   buf[11] = st_blksize, buf[12] = st_blocks
 *
 * Returns 0 on success, errno value on failure.
 */
int32_t stat_(const char *name, int32_t *buf, int64_t name_len)
{
    struct stat st;
    char *fname;
    int rc;

    if (name_len <= 0)
        return -1;

    fname = (char *)malloc((size_t)name_len + 1);
    if (!fname)
        return -1;

    memcpy(fname, name, (size_t)name_len);
    fname[name_len] = '\0';

    rc = stat(fname, &st);
    free(fname);

    if (rc != 0)
        return rc;

    buf[0]  = (int32_t)st.st_dev;
    buf[1]  = (int32_t)st.st_ino;
    buf[2]  = (int32_t)st.st_mode;
    buf[3]  = (int32_t)st.st_nlink;
    buf[4]  = (int32_t)st.st_uid;
    buf[5]  = (int32_t)st.st_gid;
    buf[6]  = (int32_t)st.st_rdev;
    buf[7]  = (int32_t)st.st_size;
    buf[8]  = (int32_t)st.st_atime;
    buf[9]  = (int32_t)st.st_mtime;
    buf[10] = (int32_t)st.st_ctime;
    buf[11] = (int32_t)st.st_blksize;
    buf[12] = (int32_t)st.st_blocks;

    return 0;
}
