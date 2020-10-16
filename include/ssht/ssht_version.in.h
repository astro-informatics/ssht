#ifndef SSHT_VERSION_H
#define SSHT_VERSION_H
inline const char *ssht_version_string() { return "@PROJECT_VERSION@"; }
inline const char *ssht_info() {
  return "package:\n"
         "  name: SSHT\n"
         "  description: Fast and accurate spin spherical harmonic transforms\n"
         "  authors:\n"
         "      - Jason McEwen\n"
         "      - Chris Wallis\n"
         "      - Martin Buttner\n"
         "      - Boris Leistedt\n"
         "      - Yves Wiaux\n"
         "  license: GPL-3\n"
         "  url: https://astro-informatics.github.io/ssht\n"
         "  version: @PROJECT_VERSION@\n";
};
// clang-format off
inline int ssht_version_major() { return @PROJECT_VERSION_MAJOR@; }
inline int ssht_version_minor() { return @PROJECT_VERSION_MINOR@; }
inline int ssht_version_patch() { return @PROJECT_VERSION_PATCH@; }
#define SSHT_VERSION_MAJOR  @PROJECT_VERSION_MAJOR@
#define SSHT_VERSION_MINOR  @PROJECT_VERSION_MINOR@
#define SSHT_VERSION_PATCH  @PROJECT_VERSION_PATCH@
// clang-format on
#endif
