Submission of ppgam_1.0.2

********

This update mainly includes a bug fix to remove a set.seed() call, which was meant just for testing, and some more S3 methods added.

********

## R CMD check results

0 errors | 0 warnings | 1 note

The following NOTE comes up, but I'm sure it's a valid DOI.

* checking CRAN incoming feasibility ... NOTE

Found the following (possibly) invalid DOIs:
  DOI: 10.1080/01621459.2016.1180986
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

********

Checked on the following:
* R 4.2.1 on Ubuntu 20.04 with gcc compiler
and also checked with rhub v2 on the following, which we all fine except 'rchk' with the error 'No files were found with the provided path: check. No artifacts will be uploaded.':

 1 [VM] linux          R-* (any version)                     ubuntu-latest on GitHub
 2 [VM] macos          R-* (any version)                     macos-13 on GitHub
 3 [VM] macos-arm64    R-* (any version)                     macos-latest on GitHub
 4 [VM] windows        R-* (any version)                     windows-latest on GitHub
 5 [CT] atlas          R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)
 6 [CT] c23            R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
 7 [CT] clang-asan     R-devel (2024-07-22 r86913)           Ubuntu 22.04.4 LTS
 8 [CT] clang16        R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
 9 [CT] clang17        R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
10 [CT] clang18        R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
11 [CT] clang19        R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
12 [CT] donttest       R-devel (2024-07-20 r86909)           Ubuntu 22.04.4 LTS
13 [CT] gcc13          R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)
14 [CT] gcc14          R-devel (2024-07-22 r86913)           Fedora Linux 40 (Container Image)
15 [CT] intel          R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)
16 [CT] mkl            R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)
17 [CT] nold           R-devel (2024-07-22 r86913)           Ubuntu 22.04.4 LTS
18 [CT] nosuggests     R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)
19 [CT] rchk           R-devel (2024-07-22 r86913)           Ubuntu 22.04.4 LTS
20 [CT] ubuntu-clang   R-devel (2024-07-22 r86913)           Ubuntu 22.04.4 LTS
21 [CT] ubuntu-gcc12   R-devel (2024-07-22 r86913)           Ubuntu 22.04.4 LTS
22 [CT] ubuntu-next    R-4.4.1 (patched) (2024-07-08 r86913) Ubuntu 22.04.4 LTS
23 [CT] ubuntu-release R-4.4.1 (2024-06-14)                  Ubuntu 22.04.4 LTS
24 [CT] valgrind       R-devel (2024-07-22 r86913)           Fedora Linux 38 (Container Image)

********

Best wishes,

Ben Youngman.


