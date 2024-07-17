A  common challenge in multimedia data understanding is the unsupervised discovery of recurring patterns, or motifs, in time series data. 

The discovery of motifs in uni-variate time series is a well studied problem and, while being a relatively new area of research, 
there are also several proposals for multi-variate motif discovery. Unfortunately, motif search among multiple variates is an expensive process, 
as the potential number of sub-spaces in which a pattern can occur increases exponentially with the number of variates. 

Consequently, many multi-variate motif search algorithms make simplifying assumptions, such as searching for motifs across all variates individually, 
assuming that the motifs are of the same length, or that they occur on a fixed subset of variates. In this paper, we are interested in addressing a relatively broad form of multi-variate motif detection, 
which seeks frequently occurring patterns (of possibly differing lengths) in sub-spaces of a multi-variate time series. 

The basic Idea behind SMM is that leveraging contextual information to help select contextually salient patterns and identify the most frequent patterns among all. 

Based on these goals, we first define the contextually salient multi-variate motif (𝙲𝚂-𝚖𝚘𝚝𝚒𝚏) discovery problem and then 

Realize a salient multi-variate motif (SMM) algorithm that, unlike existing methods, is able to seek a broad range of patterns in multi-variate time series.

Definition 11 (𝙲𝚂
-𝚖𝚘𝚝𝚒𝚏
). Let the set 𝑟𝑒𝑐𝑃𝑎𝑡={𝕊0,…,𝕊𝑘−1}
 be a (𝑘,𝑙,𝜖)
-recurring CS-pattern of the multi-variate time series 𝕋
. We refer to this as a(𝑘,𝑙,𝜖)
 𝙲𝚂
-𝚖𝚘𝚝𝚒𝚏
if the set is maximal, i.e., there is no (𝑘′,𝑙,𝜖)
-recurring CS-pattern 𝑟𝑒𝑐𝑃𝑎𝑡′
 such that 𝑘′>𝑘
 and 𝑟𝑒𝑐𝑃𝑎𝑡′⊃𝑟𝑒𝑐𝑃𝑎𝑡
.
Note that, depending on how the distance function Δ𝐶
 deals with variate alignment (e.g., whether it enforces that the sub-sequences have the same variate scopes and how it takes into account the contextual information C), we can obtain different types of 𝙲𝚂
-𝚖𝚘𝚝𝚒𝚏
s 
