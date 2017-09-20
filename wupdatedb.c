#include <ftw.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "sada_csa/interface.h"
#include "buf.h"
#include "die.h"

#define die_index(msg, err)             die("%s: %s", msg, error_index(err))

static struct {
  void   *sa;
  char   *buf;
  size_t *stack;
  size_t *offset;
  size_t *parent;
} g;

static int process_entry(const char *path, const struct stat *st, int type, struct FTW *ftw) {
  const char *p, *basename;

  basename = path + ftw->base;
  if (ftw->level == 0) {
    basename = path;
    buf_push(g.buf, '\n');
  }

  if (buf_len(g.stack) > ftw->level)
    buf_resize(g.stack, ftw->level);

  if (type == FTW_D)
    buf_push(g.stack, buf_len(g.offset));

  buf_push(g.offset, buf_len(g.buf));
  buf_push(g.parent, ftw->level > 0 ? g.stack[ftw->level - 1] : 0);

  for (p = basename; *p; p++)
    buf_push(g.buf, *p);
  buf_push(g.buf, '\n');

  return 0;
}

static int process_dir(const char *path) {
  static const int NOPENFD = 512;
  int ret;

  ret = nftw(path, process_entry, NOPENFD, 0);
  buf_push(g.buf, '\0');
  return ret;
}

static int write_opi(const char *filename) {
  FILE *f;
  char path[PATH_MAX];

  snprintf(path, sizeof(path), "%s.opi", filename);
  f = fopen(path, "w");
  if (!f)
    return -1;

  fwrite(g.offset, sizeof(*g.offset), buf_len(g.offset), f);
  fwrite(g.parent, sizeof(*g.parent), buf_len(g.parent), f);
  fclose(f);

  return 0;
}

static void display_usage(void) {
  fprintf(stderr, "usage: wupdatedb [OPTIONS] root\n");
  fprintf(stderr, "  options:\n");
  fprintf(stderr, "    -o output filename (You need to exclude the extension)\n");
  fprintf(stderr, "    -h displays this message\n");
}

int main(int argc, char *argv[]) {
  const char *root = ".";
  const char *filename = NULL;
  char build_options[PATH_MAX];
  int c, err;

  opterr = 0;
  while ((c = getopt(argc, argv, "o:h?")) != -1) {
    switch (c) {
      case 'o':
        filename = optarg;
        break;
      case 'h':
      case '?':
        display_usage();
        return 0;
    }
  }

  if (filename == NULL) {
    display_usage();
    return -1;
  }

  if (optind < argc)
    root = argv[argc - 1];

  if (process_dir(root))
    die_errno("process_dir");

  if (write_opi(filename))
    die_errno("write_opi");

  snprintf(build_options, sizeof(build_options), "filename=%s", filename);
  if ((err = build_index((unsigned char *)g.buf, buf_len(g.buf), build_options, &g.sa)))
    die_index("build_index", err);

  if ((err = free_index(g.sa)))
    die_index("free_index", err);

  buf_free(g.buf);
  buf_free(g.stack);
  buf_free(g.offset);
  buf_free(g.parent);

  return 0;
}
