# Null covariate generation manifest
Generated: 2026-09-16 10:14:05

Omega source: `module_manuscript_rho05/baypass_stage1/aland_excluded/omega_mat_omega.out`
Omega MD5: `7abb549341c5c08a7ba238892c7c0740`
Omega max asymmetry (pre-symmetrizing): 0.000e+00
Omega negative eigenvalues clamped: 0 / 19

Population order (P=19, source: data/hybrids_only_maf005.Rdata, Aland excluded):
  Bunkkeri, Tvarminne, Grundsund, Heinamaki, Hiivola, Jarvenpaa, Karsikas, Katiskoski, Kummunmaki, LangholmenR, LangholmenW, Nyrhispera74, Nyrhispera75, Parikkala, Pikkala, Sielva, Svanvik1, Svanvik2, Vuosaari

Real mitoC2 split (source: u.mito_contrast): 7 = +1 (Faquilonia-like), 12 = -1 (Fpolyctena-like)

Replicates: 10
Structured seeds: 20101, 20102, 20103, 20104, 20105, 20106, 20107, 20108, 20109, 20110
Permute seeds: 20201, 20202, 20203, 20204, 20205, 20206, 20207, 20208, 20209, 20210

All validation checks in R/01_generate_null_covariates.R passed (mean/sd of
structured draws; permutation integrity of unstructured draws; exact
7/12 group sizes for both mitoC2 null variants).
