# Function matrix of appartenance group

Method to create a binary matrix with proteins in columns and peptides
in lines on a `MSnSet` object (peptides)

## Usage

``` r
BuildAdjacencyMatrix(obj.pep)
```

## Arguments

- obj.pep:

  An object (peptides) of class `MSnSet`.

## Value

A binary matrix

## Author

Florence Combes, Samuel Wieczorek, Alexia Dorffer

## Examples

``` r
data(subR25pept)
BuildAdjacencyMatrix(subR25pept[[1]])
#> 100 x 96 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 96 column names ‘1005’, ‘1017’, ‘103’ ... ]]
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . 1 . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . 1 1 . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . 1 . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . 1 . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . 1 . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . 1 . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . 1
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . 1 . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . 1 . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . 1 . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . 1 . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . 1 . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . 1 . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          1 . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . 1 . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . 1 . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . 1 .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . 1 . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . 1 . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . 1 . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . 1 . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . 1 . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . 1 . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               1 . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . 1 . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . 1 . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . 1 . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . 1 . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . 1
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . 1 .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . 1 . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . 1 . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . 1 . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . 1 . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . 1 . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . 1 . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . 1 . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . 1 . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . 1 . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    1 . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . 1 . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . 1 .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . 1 . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . 1 . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . 1 . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . 1 . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . 1 . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . 1 . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . 1 . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . 1 . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . 1 . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . 1 . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . 1 . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . 1
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . 1 . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . 1 . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . 1 . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . 1 . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . 1 . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . 1 . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . 1 . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . 1 . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . 1 . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . 1 . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . 1
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . 1 . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . 1 . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . 1 . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . 1 . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . 1 . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . 1 . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . 1 .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . 1 . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . 1 . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . 1 . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . 1 . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               1 . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . 1 . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . 1 . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . 1 . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . 1 . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . 1 . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . 1 .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . 1 . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . 1 . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . 1 . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . 1 . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . 1 . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . 1 . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . 1 . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . 1 . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . 1
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . 1 . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . 1 . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . 1 . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . 1 . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . 1 . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          1 . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . 1 . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . 1 . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . 1 . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                           
#> AAAAQDEITGDGTTTVVCIVGEIIR                .
#> AAADAISDIEIK                             .
#> AAADAISDIEIKDSK                          .
#> AAAEEQAKR                                .
#> AAAEGVANIHIDEATGEMVSK                    .
#> AAAEYEKGEYETAISTINDAVEQGR                .
#> AAAHSSIKEYDQAVK                          .
#> AAAPGIQIVAGEGFQSPIEDR                    .
#> AAAPTVVFIDEIDSIAK                        .
#> AACIVQNGIATWFPIAVTK                      .
#> AADAIIIK                                 .
#> AADETAAAFYPSK                            .
#> AADIINIAK                                .
#> AADIPVVGNAAGHSNDWFDIK                    .
#> AADTPETSDAVHTEQKPEEEKETIQEE              .
#> AAEAATTDITYR                             .
#> AAEAGETGAATSATEGDNNNNTAAGDK              .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             .
#> AAEEADADAEIADEEDAIHDEI                   .
#> AAEIDVINDPK                              .
#> AAEIIIENR                                .
#> AAEIIISDQDNVIPK                          .
#> AAENASNAIAETR                            .
#> AAESIISIANVPDGDSR                        .
#> AAFDEDGNISNVK                            .
#> AAFISAIVGK                               .
#> AAFNGVTFK                                .
#> AAFNYQFDSIIEHSEK                         .
#> AAFTECCQAADK                             .
#> AAFTYIINDPEIAK                           .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK .
#> AAGANVDNVWADVYAK                         .
#> AAGFIIEK                                 .
#> AAGGVVEIIA                               .
#> AAGIFVSTSSYGGGQESTVK                     .
#> AAGITAAYAR                               .
#> AAGIVDIATVISTSAYIIESK                    .
#> AAGKEIGDFEDISTENEK                       .
#> AAIANVYEYR                               .
#> AAIDFYTK                                 .
#> AAIEAGAFEAVTSNHWAEGGK                    .
#> AAIEDGPVTAENISSETAR                      .
#> AAIEDGWVPGK                              .
#> AAIEEIVK                                 .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            .
#> AAIGSSPINFPSSSQR                         .
#> AAIIACAAEYIQK                            .
#> AAIIGSIGSIFK                             .
#> AAIINQYFAQAYK                            .
#> AAIISSGNVK                               .
#> AAINDPAKAPIIINNIIDSGIR                   .
#> AAIQTYIPK                                .
#> AAISFGAKPEEQK                            .
#> AAITDFER                                 .
#> AAITIIQFDGTGTR                           .
#> AAIVQIDATPFR                             .
#> AAIYAIHSIGCK                             .
#> AANAIKDIYGWTQTSIDDYPIK                   .
#> AANAPVYVTR                               .
#> AANIAHDNQTTVEAYK                         .
#> AANIGGVAVSGIEMAQNSQR                     .
#> AANQGAIPPDISIIVK                         .
#> AANQTASSIVDFYNAIGDDEEEK                  .
#> AANSHRIIDIQESQANCSHFFIEPIK               .
#> AAPIVDDEETEFDIYNSK                       .
#> AAPSPISHVAEIR                            .
#> AAQDVWNR                                 .
#> AAQIGFNTACVEK                            .
#> AAQIGSSFIAQIK                            .
#> AAQQQWGNDFYK                             .
#> AASDAIPPASPK                             .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              .
#> AASKPFIETFICGR                           .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          .
#> AASSINRVDTIR                             .
#> AASYADKINTPEIWSQIGTAQIDGIR               .
#> AATIIPQFVGIK                             .
#> AATIISNII                                .
#> AATVVATSDCIIWAIDR                        .
#> AAVDCECEFQNIEHNEK                        .
#> AAVEEGIIPGGGTAIVK                        .
#> AAVPFNREQIESVIR                          .
#> AAVSGKPYFFFGSDSAPHPVQNK                  .
#> AAWWSPTGDYIAFIK                          .
#> AAYAIGGIGSGFANNEK                        .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    .
#> AAYGDISDEEEK                             .
#> AAYSYMFDSIR                              .
#> ACAAQTNATFIK                             .
#> ACDTSNDNFPIQYDGSK                        .
#> ACGIFSGYPDTFK                            1
#> ACGIIISEER                               .
#> ACGVSRPVIAASITTNDASAIK                   .
#> ACNFQFPEIAYPGK                           .
#> ACPVGNEAGVTTSIR                          .
#> ACVIVVSDIK                               .
#> ACVVYGGSPIGNQIR                          .
#> ADAEWVQSTASK                             .
#> ADASGEGVEDEASGVHK                        .
```
