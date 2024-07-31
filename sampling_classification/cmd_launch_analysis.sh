### build datasets

module load r/4.1.1
sbatch --mem-per-cpu=30GB -A myProject -t 12:00:00 --wrap="Rscript --vanilla build_datasets.R"

### launch classification for fold 1 sampling strategy hBA, level subphylum, 
### classification method RandomForest and no computation of shapley values

k=1 # nb fold: 1 to 4
s=smote #  V1 / V3 / smote
l=subphylum # subphylum class order superfamily genus species
m=RF # svm rapide knn RF
sh=FALSE # TRUE ou FALSE

sbatch --mem-per-cpu=40GB -A myProject -c 6 --export=fold=k,sampling=s,level=l,method=m,shap=sh run_classif.sh

