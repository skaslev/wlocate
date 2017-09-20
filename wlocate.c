#include <fcntl.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "sada_csa/interface.h"
#include "sada_csa/comparray4.h"
#include "buf.h"
#include "die.h"

#define die_index(msg, err)             die("%s: %s", msg, error_index(err))

static struct {
  CSA    *sa;
  size_t  size;
  size_t  len;
  size_t *offset;
  size_t *parent;
} g;

static int read_opi(const char *filename) {
  int fd, ret = -1;
  struct stat st;
  void *m;
  char path[PATH_MAX];

  snprintf(path, sizeof(path), "%s.opi", filename);
  fd = open(path, O_RDONLY);
  if (fd < 0)
    return fd;

  if (fstat(fd, &st))
    goto out;

  m = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (!m)
    goto out;

  g.size = st.st_size;
  g.len = st.st_size / (2 * sizeof(*g.offset));
  g.offset = m;
  g.parent = g.offset + g.len;
  ret = 0;

out:
  close(fd);
  return ret;
}

static char *get_full_path(char *buf, size_t node) {
  long ch, pos;
  const size_t parent = g.parent[node];

  if (parent != node) {
    buf = get_full_path(buf, parent);
    if (buf_len(buf) && buf_last(buf) != '/')
      buf_push(buf, '/');
  }

  for (pos = csa_inverse(g.sa, g.offset[node] + 1);
       (ch = csa_T(g.sa, pos)) != '\n';
       pos = csa_psi(g.sa, pos))
    buf_push(buf, ch);

  return buf;
}

static int size_cmp(const void *lhs, const void *rhs) {
  const size_t l = *((const size_t *)lhs);
  const size_t r = *((const size_t *)rhs);
  return l - r;
}

static void do_search(const char *pattern) {
  static const unsigned long PADDING = 256;
  int err;
  char *buf = NULL;
  unsigned char *snippet_text;
  unsigned long i, pattern_len, match_len;
  unsigned long numocc, *snippet_len, *abs_from, *abs_occ;

  pattern_len = strlen(pattern);
  match_len = pattern_len + 2 * PADDING;
  if ((err = display(g.sa, (const unsigned char *)pattern, pattern_len,
                     PADDING, &numocc, &snippet_text, &snippet_len,
                     &abs_from, &abs_occ)))
    die_index("display", err);

  for (i = 0; i < numocc; i++) {
    size_t key, *match;
    unsigned char *snippet = snippet_text + match_len * i;
    unsigned char *start = snippet + (abs_occ[i] - abs_from[i]);

    while (start > snippet && *start != '\n')
      start--;
    if (*start == '\n')
      start++;

    key = abs_from[i] + (start - snippet);
    match = bsearch(&key, g.offset, g.len, sizeof(*g.offset), size_cmp);
    if (!match) {
      fprintf(stderr, "WTF failed to find offset %zu\n", key);
      continue;
    }

    buf_resize(buf, 0);
    buf = get_full_path(buf, match - g.offset);
    buf_push(buf, '\n');
    fwrite(buf, 1, buf_len(buf), stdout);
  }

  buf_free(buf);
  if (numocc) {
    free(abs_from);
    free(abs_occ);
    free(snippet_len);
    free(snippet_text);
  }
}

static void display_usage(void) {
  fprintf(stderr, "usage: wlocate [OPTIONS] pattern\n");
  fprintf(stderr, "  options:\n");
  fprintf(stderr, "    -d database filename to use (You need to exclude the extension)\n");
  fprintf(stderr, "    -h displays this message\n");
}

int main(int argc, char *argv[]) {
  const char *filename = NULL;
  const char *pattern = NULL;
  int c, err;

  opterr = 0;
  while ((c = getopt(argc, argv, "d:h?")) != -1) {
    switch (c) {
      case 'd':
        filename = optarg;
        break;
      case 'h':
      case '?':
        display_usage();
        return 0;
    }
  }

  if (!filename || optind != argc - 1) {
    display_usage();
    return -1;
  }

  pattern = argv[argc - 1];

  if (read_opi(filename))
    die_errno("read_opi");

  if ((err = load_index(filename, (void**)&g.sa)))
    die_index("load_index", err);

  do_search(pattern);

  if ((err = free_index(g.sa)))
    die_index("free_index", err);

  if (munmap(g.offset, g.size))
    die_errno("munmap");

  return 0;
}
