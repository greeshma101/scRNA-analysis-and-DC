# 👋 Hi, I'm Greeshma Gopinath  

🎓 Aspiring PhD researcher in **bioinformatics and neurodegeneration**  
🔬 Researching **cell-to-cell communication in glaucoma** through **scRNA-seq & proteomics**  
🌱 Passionate about integrating **statistics, omics data, and biology** to uncover molecular mechanisms of disease  

---

## 📘 My Research Journey  

I began my scientific journey with a strong **statistical background** during my Master’s degree.  
Through coursework in **RStudio**, I learned methods such as:  
- Principal Component Analysis (PCA)  
- Mixed linear models  
- Regression approaches  

This gave me a quantitative foundation and became a **stepping stone into bioinformatics**.  

As biology has become increasingly **data-driven and omics-focused**, I grew fascinated by the idea that understanding molecular mechanisms begins at the **cellular level**. That curiosity drove me to teach myself new computational methods and dive deeper into large-scale biological data analysis.  

---

## 🔬 Current Research Focus  

My research centers on understanding **how retinal cells communicate and degenerate in glaucoma**.  



### 🧬 What I’ve Explored

- **Single-cell RNA-seq analysis (scRNA-seq)** of the retina and anterior eye  
  I primarily **re-analyzed public scRNA-seq datasets**, which was a fantastic way to **gain hands-on experience and deepen my understanding** of the analysis workflow. These datasets also allowed me to **validate my findings**. In particular, I worked with the dataset from **Fadl et al.2020**  

<p align="center">
  <img src="figures/umap_retina.png" alt="UMAP of scRNA retina data" width="45%"/>
  <img src="figures/deconv_ev.png" alt="Deconvolution of EV proteomics data" width="45%"/>
</p>

- **Proteomics-based deconvolution**  
  - Built a **cellular reference matrix** from scRNA-seq cell types  
  - Traditional deconvolution methods are usually applied to **bulk transcriptomics data**, but I adapted a **proteomics-based approach** inspired by [Newman et al., 2015](https://doi.org/10.1038/nmeth.3337) using the **SVR method**  
  - Applied this to **extracellular vesicle (EV) proteomics data** to estimate contributions of different retinal cell populations  
  - This strategy addresses the **challenge of isolating EVs by cell type experimentally**, providing a computational solution to dissect cellular contributions
