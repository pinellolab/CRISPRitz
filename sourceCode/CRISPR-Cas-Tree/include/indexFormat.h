#pragma once
// CRISPRitz TST index (.bin) format constants. The layout, and what a format-2
// index holds, are specified in docs/INDEX_FORMAT.md.
//
// Included by path ("include/indexFormat.h") from both mainParallel.cpp (buildTST)
// and searchOnTST.cpp (searchTST): the buildTST recipe passes no -I for this
// directory.

// A format-2 .bin starts, and ends, with these 4 bytes. Read as a little-endian
// int32 they are -28747197: a format-1 reader takes the first int as the leaf
// count, and refuses a non-positive one (CRISPRitz 2.9.1) or fails to allocate it
// (earlier), so an older engine fails loudly on a format-2 index instead of
// misreading it.
static const char CZI_MAGIC[4] = {'C', 'Z', 'I', (char)0xFE};
static const int CZI_FORMAT_VERSION = 2;
// magic + version + header_len + flags + guide_len + pam_len + max_bulges, then
// the pam_len characters of the build PAM region
static const int CZI_HEADER_FIXED = 28;

// header flags; a reader refuses a file with a bit it does not know
static const int CZI_FLAG_PAM_AT_START = 1 << 0; // the PAM is at the 5' end
static const int CZI_FLAG_PAMLESS = 1 << 1;      // the build PAM is all N
static const int CZI_FLAG_FOUR_ALLELE = 1 << 2;  // the tree and the flanks hold Q (below)
static const int CZI_FLAG_KNOWN = CZI_FLAG_PAM_AT_START | CZI_FLAG_PAMLESS | CZI_FLAG_FOUR_ALLELE;

// Q, the fifth symbol: a position whose reference base and SNP ALTs cover A, C, G
// and T (written by the enricher), matched as any base (1111), where N stays a gap
// (0000) that matches nothing. Each 4-bit slot of the tree already uses 15 of its 16
// values -- a root, lo or hi slot: 0 = no kid and the 14 IUPAC letters; an eq slot:
// the 14 letters and 15 = end of string -- so Q takes the value its slot never uses:
// 15 in a root, lo or hi slot (its A|C|G|T mask), 0 in an eq slot. A flank nibble Q
// is 15. Only a file flagged CZI_FLAG_FOUR_ALLELE holds them; a file without Q is
// written byte for byte as before, so an older reader refuses only what it cannot
// decode.
static const unsigned char CZI_NIBBLE_Q = 0x0F;    // root, lo and hi slots; flanks
static const unsigned char CZI_NIBBLE_Q_EQ = 0x00; // eq slots

// per-leaf flags byte, which takes the place of format 1's '0'/'_' marker
static const unsigned char CZI_LEAF_HAS_NEXT = 0x80;  // an int32 next follows
static const unsigned char CZI_LEAF_PAM_LEFT = 0x40;  // window = flank|guide|pad
static const unsigned char CZI_LEAF_FLANK_OFF = 0x20; // part of the flank is off the contig
static const unsigned char CZI_LEAF_PAD_MASK = 0x1F;  // pad_valid, 0..31

// pad_valid has 5 bits, so an index cannot be built with more bulges than this
static const int CZI_MAX_BULGES = 31;
