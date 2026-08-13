#
#
#
#
# UGH. Hahah


library("data.table")

effect_formulas <- fread("../metatools_dev/data/effect_size_formulas.csv")
effect_formulas

effect_formulas[escalc_name == "", ]

unique(effect_formulas$effect_size)
effect_formulas <- effect_formulas[effect_size %in% c("lnRoM", "SMD", "lnOR", "lnRR", "lnCVR",
                                                       "lnHWE_A", "lnVR", "SMDH", "Zr", "lnM")]

effect_formulas[effect_size == "lnHWE_A", effect_size := "lnHWE"]

names(effect_formulas)
effect_formulas <- effect_formulas[, !c("lower_filter", "upper_filter", "n_versions")]

# file.exists("../metatools/data/effect_formulas.rda")
# save(effect_formulas, file = "../metatools/data/effect_formulas.rda")

# Conversion formulas:

conversion_formulas <- fread("../metatools_dev/data/conversion_formulas.csv")
conversion_formulas

# save(conversion_formulas, file = "../metatools/data/conversion_formulas.rda")

usethis::use_data(effect_formulas, conversion_formulas, internal = FALSE, overwrite = TRUE)
