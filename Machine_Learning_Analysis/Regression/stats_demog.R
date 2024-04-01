

HC=filter(data_matched, Diagnosis==0)
SCZ=filter(data_matched, Diagnosis==1)
HC=subset(HC, SITE_ID %in% c(1,8,9))
SCZ=subset(SCZ, SITE_ID %in% c(1,8,9))
dim(HC)
dim(SCZ)

mean(HC$AGE)
mean(SCZ$AGE)
sd(HC$AGE)
sd(SCZ$AGE)
t_test_result <- t.test(HC$AGE, SCZ$AGE)

# Extract p-value
p_value <- t_test_result$p.value
p_value

sex_female=filter(HC,SEX==1)
dim(sex_female)

sex_female=filter(SCZ,SEX==1)
dim(sex_female)


# Assuming HC and SZ are data frames containing the relevant data

# Frequency table for males in HC and SZ groups
male_HC <- sum(HC$SEX == 0)
male_SZ <- sum(SCZ$SEX == 0)

# Frequency table for females in HC and SZ groups
female_HC <- sum(HC$SEX == 1)
female_SZ <- sum(SCZ$SEX == 1)

# Combine frequency tables into a contingency table
contingency_table <- matrix(c(male_HC, male_SZ, female_HC, female_SZ), nrow = 2, byrow = TRUE)
rownames(contingency_table) <- c("Male", "Female")
colnames(contingency_table) <- c("HC", "SZ")

# Print contingency table
print(contingency_table)

# Calculate chi-square test
chi_square_result <- chisq.test(contingency_table)

# Print chi-square test result
print(chi_square_result)
# Assuming HC and SZ are data frames containing the relevant data

# Create contingency table

contingency_table <- table(Group = c(rep("HC", sum(HC$SEX == 0)), rep("SCZ", sum(SCZ$SEX == 0), 
                                                                      Gender = c(rep("Male", sum(HC$SEX == 0)), rep("Male", sum(SZ$SEX == 0)))
)

# Calculate chi-square test
chi_square_result <- chisq.test(contingency_table)

# Print chi-square test result
print(chi_square_result)


                                                                      

