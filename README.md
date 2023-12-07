# PCC_Mutations

Attached in this repository are the codes and data used in our publication *Uncovering potential interventions for pancreatic cancer patients via mathematical modeling*[^1][^2]. Specifically, you will find:

TCGA pancreatic cancer data: 
- expressions -- gene expression data
- mutation_data_public -- aggregate statistics of expression data, along with mutation attractor information
- name_map -- a mapping of gene names with our their corresponding nodes in our model
- paad_mutation_status_2018july26 -- time-course data
- PCC_expressionsRed -- subset of the original data focusing on four main mutation combinations

<br />

Statistical analysis code (Python3)
- KaplanMeier_PCC -- survival curves
- PC_KW_comb -- Kruskal-Wallis test and plotting 

<br />

Model functions
- pcc_logical_functions -- functions used in the model
  - use these functions for stable motifs (https://github.com/jgtz/StableMotifs)
    - can be used to derive attractors and controls

<br />

Simulation code (MATLAB)
- comb_pcc_DP -- cumulative code to run simulations using SDDS
    - instructions commented throughout code
- MatlabToolboxes -- folder of all the needed toolboxes to run MATLAB code

<br />

Computational algebra control file
- pcc_mutation_node_control_public -- use with (https://www.unimelb-macaulay2.cloud.edu.au/)
- pcc_mutation_edge_control -- use with (https://www.unimelb-macaulay2.cloud.edu.au/)
  - instructions commented throughout code

<br />

We have included a *How-To* file that details how to find and simulate controls
- How-To

<br />

NOTICE: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


[^1]:  Uncovering potential interventions for pancreatic cancer patients via mathematical modeling. Daniel Plaugher, Boris Aguilar, David Murrugarra. *Journal of Theoretical Biology*, accepted, 2022. https://authors.elsevier.com/c/1fJ4457imC48X
[^2]: Coming soon!
