// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Constants used for PreMeta.

#ifndef PREMETA_CONSTANTS_H
#define PREMETA_CONSTANTS_H

#include <stdint.h>

namespace premeta {

// General Constants.
extern const int64_t kMaxBufferSize;
extern const char kAlleleFilePositionHeader[];
extern const char kAlleleFileSnpIdHeader[];
extern const char kAlleleFileMetaSnpIdHeader[];
extern const char kAlleleFileStudySnpIdHeader[];
extern const char kAlleleFileMajorHeader[];
extern const char kAlleleFileMinorHeader[];

// R Constants.
extern const int kRXDRIntegerBytes;
extern const int kRXDRDoubleBytes;
extern const int kItemEndMarker;

// Allele File Constants.
extern const char kAlleleFileDefaultDelimiter[];
extern const char kAlleleFileDefaultCommentChar[];

// RareMetal Constants.
extern const char kRareMetalScoreFileDefaultDelimiter[];
extern const char kRareMetalScoreFileDefaultCommentChar[];
extern const char kRareMetalScoreFileSampleNum[];
extern const char kRareMetalGroupFileDefaultDelimiter[];
extern const char kRareMetalGroupFileDefaultCommentChar[];
extern const char kRareMetalCovFileDefaultDelimiter[];
extern const char kRareMetalCovFileDefaultCommentChar[];
extern const char kResidualVarianceStr[];
extern const char kDummyRefAllele[];
extern const char kDummyAltAllele[];
extern const int kRareMetalWindowSize;

// Mass Constants.
extern const char kMassDefaultDelimiter[];
extern const char kMassDefaultCommentChar[];
extern const char kSigmaSquared[];

// MetaSKAT Constants.
extern const char kMInfoDefaultDelimiter[];
extern const char kMInfoDefaultCommentChar[];
extern const char kMInfoHeaderSetId[];
extern const char kMInfoHeaderSetIdNumeric[];
extern const char kMInfoHeaderSnpId[];
extern const char kMInfoHeaderScore[];
extern const char kMInfoHeaderMaf[];
extern const char kMInfoHeaderMissingRate[];
extern const char kMInfoHeaderMajorAllele[];
extern const char kMInfoHeaderMinorAllele[];
extern const char kMInfoHeaderPass[];
extern const char kMInfoHeaderStartPos[];
extern const int kMaxMatrixSize;
extern const uint32_t crc32_table[];

}  // namespace premeta
#endif
