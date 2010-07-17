#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>

int main (int argc, char *argv[]) {
    int status;
    struct stat st_buf;

    // Ensure argument passed.

    if (argc != 2) {
        printf ("Usage: progName <fileSpec>\n");
        printf ("       where <fileSpec> is the file to check.\n");
        return 1;
    }

    // Get the status of the file system object.

    status = stat (argv[1], &st_buf);
    if (status != 0) {
        printf ("Error, errno = %d\n", errno);
        return 1;
    }

    // Tell us what it is then exit.

    if (S_ISREG (st_buf.st_mode)) {
        printf ("%s is a regular file.\n", argv[1]);
    }
    if (S_ISDIR (st_buf.st_mode)) {
        printf ("%s is a directory.\n", argv[1]);
    }

    return 0;
