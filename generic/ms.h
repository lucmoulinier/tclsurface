
typedef struct {
  int   Id;              /* Id used to restore F<Dim>F                   */
  int   LengthOfHeader;  /* Length of header in bytes before any vector  */
  int   Dim;             /* Number of dimensions                         */
  int   NoOfSamples;     /* Total number of samples                      */
  int   Type;		 /* Picture type (eg. format on disc, domain)	 */
  char  FileName[160];	 /* Name used for saving/restoring this image	 */
  char  Description[80]; /* eg. short description of image (optional)	 */
  char  Date[8];	 /* Date (DDMMYY) (optional)			 */
  int   ArrayType;	 /* Defines number format: Complex(2) or Real(1) */
  int   *Samples;        /* Pointer to number of samples. Length is Dim  */
  float *Min;            /* Pointer to min positions. Length is Dim      */
  float *Delta;          /* Pointer to sampling distances. Length is Dim */
  float *Signal;         /* Pointer to signal. Length: product of dims   */
} MultiStruct;								  
