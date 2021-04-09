package org.bwh.pathology;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SampleIdentityChecker {

    public static void main(String[] args) {

        Map<String, Map<String, String>> sample_genotypes = new HashMap<String, Map<String, String>>();

        try {

            BufferedReader br = new BufferedReader(new FileReader(args[0]));

            String line;

            // First load the contents of file into HashMap
            while ((line = br.readLine()) != null) {
                String[] arr = line.split("\t");
                int len = arr.length;
                //Skip  header
                if (!arr[0].equalsIgnoreCase("sample")) {
                    String sample_name = arr[0];
                    Map<String, String> geno = new HashMap<String, String>();

                    for (int i = 1; i < len; i++) {

                        String site = "site" + i;

                        geno.put(site, arr[i]);
                    }
                    sample_genotypes.put(sample_name, geno);

                }

            }

            //Hash set to keep track of samples alreadyh compared to avoid duplication
            Set<String> samples_processed = new HashSet<String>();


            for (String sample : sample_genotypes.keySet()) {
                for (String testsample : sample_genotypes.keySet()) {
                    if (!sample.equalsIgnoreCase(testsample) && !samples_processed.contains(testsample)) {

                        compareGenotypes(sample, testsample, sample_genotypes.get(sample), sample_genotypes.get(testsample));
                    }

                    samples_processed.add(sample);

                }


            }
            br.close();

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
    }

    /**
     * Get total read counts given SNP site coverage
     *
     * @param s
     * @return
     */

    private static int getCoverage(String s) {

        String[] arr3 = s.split(",");

        int count = 0;
        return count + Integer.parseInt(arr3[0]) + Integer.parseInt(arr3[1]);

    }

    /**
     * Compare the genotypes between 2 samples and prints out if there are to likely match
     * <p>
     * 2 Methods
     * <p>
     * Nonstatistical method
     * <p>
     * Check if the total reads delta to obtain a match between 2 samples is below or equal to error rate.
     * For examples 2 samples with cov 275 and 280 can have maximum of 6 read errors ( based on 0.01).
     * If min reads required for match between 2 samples is <=6 then there is a possibility they are from same sample
     * <p>
     * Statistical Method
     * Compute loglikelihood if genotypes at site matched (assumes they are from same).
     * Compute loglikelihood if genotypes do not match.
     * Ratio ( LL(S1)//LL(S0)) should be significantly higher.Choose threshold 5.
     * <p>
     * If it satisfies both conditions then it is "Likely Match")
     *
     * @param sample
     * @param testsample
     * @param genotype_samp1
     * @param genotype_samp2
     * @throws IOException
     */
    private static void compareGenotypes(String sample, String testsample, Map<String, String> genotype_samp1, Map<String, String> genotype_samp2) throws IOException {

        int samp1_total_coverage = 0;
        int samp2_total_coverage = 0;
        int mismatch_coverage = 0;


        Double prob_same = 0.0;
        Double prob_diff = 0.0;
        for (String site : genotype_samp1.keySet()) {
            String samp1_geno = getGenotype(genotype_samp1.get(site));
            String samp2_geno = getGenotype(genotype_samp2.get(site));
            int samp1_coverage = getCoverage(genotype_samp1.get(site));
            int samp2_coverage = getCoverage(genotype_samp2.get(site));
            samp1_total_coverage = samp1_total_coverage + samp1_coverage;
            samp2_total_coverage = samp2_total_coverage + samp2_coverage;
            double probability_sample1 = calculateProbabilityForSample(genotype_samp1.get(site));
            double probability_sample2 = calculateProbabilityForSample(genotype_samp2.get(site));
            double probability = calculateProbability(genotype_samp1.get(site), genotype_samp2.get(site));

          // If genotypes match at the site
            if (samp1_geno.equalsIgnoreCase(samp2_geno)) {

               // Check Het or Hom
                if (samp1_geno.length() == 1) {
                    prob_same += calculateHomProbability(samp1_coverage + samp2_coverage, probability);
                } else {
                    prob_same += calculateHetProbability(samp1_coverage + samp2_coverage, probability);
                }

            }
            // If missing genotypes (no coverage) in either sample at given site
            else if (samp1_geno.equalsIgnoreCase("N") || samp2_geno.equalsIgnoreCase("N")) {

                // If Sample1 has a genotype (i.e sample2 has no coverage or genotype)
                if (samp1_coverage > 0) {

                    // Read delta
                    mismatch_coverage = mismatch_coverage + samp1_geno.length();
                    if (samp1_geno.length() > 1) {
                        prob_diff += calculateHetProbability(samp1_coverage, probability_sample1);

                    } else {
                        prob_diff += calculateHomProbability(samp1_coverage, probability_sample1);

                    }


                }

                // If Sample2 has a genotype (i.e sample1 has no coverage or genotype)
                if (samp2_coverage > 0) {

                    // Read delta
                    mismatch_coverage = mismatch_coverage + samp2_geno.length();
                    if (samp2_geno.length() > 1) {
                        prob_diff = calculateHetProbability(samp2_coverage, probability_sample2);
                    } else {
                        prob_diff += calculateHomProbability(samp2_coverage, probability_sample2);

                    }

                }
            } else { //Genotypes do not match

                // Homozygous mismatch
                if (samp1_geno.length() == samp2_geno.length()) {
                    // Read delta
                    mismatch_coverage = mismatch_coverage + 2;
                    prob_diff += calculateHomProbability(samp1_coverage, probability_sample1);
                    prob_diff += calculateHomProbability(samp2_coverage, probability_sample2);

                } else {

                    // Heterozygous mismatch
                    mismatch_coverage = mismatch_coverage + 1;
                    if (samp1_geno.length() > 1) {
                        prob_diff += calculateHetProbability(samp1_coverage, probability_sample1);
                        prob_diff += calculateHomProbability(samp2_coverage, probability_sample2);
                    } else {
                        prob_diff += calculateHomProbability(samp1_coverage, probability_sample1);
                        prob_diff += calculateHetProbability(samp2_coverage, probability_sample2);
                    }
                }

            }
        }

        int readerror_threshold = (int) Math.ceil(samp1_total_coverage * 0.01) + (int) Math.ceil(samp2_total_coverage * 0.01);
        if (mismatch_coverage <= readerror_threshold && prob_same / prob_diff > 5) {
            System.out.println("Sample1Name\tSample2Name\tNoReads_Delta\tNoReads_possibleError\tSample1_total_coverage\tsamp2_total_coverage\tFingerprintStatus\tLogLikelihood_S1\tLogLikelihood_S0\tRatio");
            System.out.println(sample + "\t" + testsample + "\t" + mismatch_coverage + "\t" + samp1_total_coverage + "\t" + samp2_total_coverage + "\t" + "Likely Match" + "\t" + prob_same + "\t" + prob_diff + "\t" + prob_same / prob_diff);
        }


    }

    /**
     * Calculate Het Site Probability given coverage and probability of allele at SNP site for single sample
     *
     * @param coverage
     * @param probability
     * @return
     */

    private static double calculateHetProbability(int coverage, double probability) {
        return (coverage * Math.log(0.5) + Math.log(probability));
    }

    /**
     * Calculate Hom Site Probability given coverage and probability of allele at SNP site for single sample
     *
     * @param coverage
     * @param probability
     * @return
     */
    private static double calculateHomProbability(int coverage, double probability) {
        return (coverage * Math.log(0.99) + Math.log(probability));
    }

    /**
     * Calculate probability of an Allele at given site for a sample given read counts
     *
     * @param s_counts
     * @return
     */

    private static double calculateProbabilityForSample(String s_counts) {

        String[] array = s_counts.split(",");
        double prob = 1.0;
        int count_a = Integer.parseInt(array[0]);
        int count_b = Integer.parseInt(array[1]);


        if (count_a > 0 && count_b > 0) {


            double maf = (double) Math.min(count_a, count_b) / (count_a + count_b);
            prob = maf * (1 - maf) * 2;

        }

        return prob;
    }

    /**
     * Calculate probability for Allele given both sample read counts
     *
     * @param s1_counts
     * @param s2_counts
     * @return
     */

    private static double calculateProbability(String s1_counts, String s2_counts) {

        String[] arr_s1 = s1_counts.split(",");
        double prob = 1.0;
        int count_a = Integer.parseInt(arr_s1[0]);
        int count_b = Integer.parseInt(arr_s1[1]);

        String[] arr_s2 = s2_counts.split(",");

        count_a = count_a + Integer.parseInt(arr_s2[0]);
        count_b = count_b + Integer.parseInt(arr_s2[1]);

        if (count_a > 0 && count_b > 0) {


            double maf = (double) Math.min(count_a, count_b) / (count_a + count_b);
            prob = maf * (1 - maf) * 2;

        }


        return prob;


    }


    /**
     * Get genotype given Allele counts from File
     *
     * @param s (0,6)
     * @return String ( A,B,AB or N)
     * <p>
     * Passing 0,6 should return B (0 Reads supporting A Allele and 6 Reads supporting B Allele)
     */
    private static String getGenotype(String s) {

        String[] array = s.split(",");

        String genotype = "";

        if (Integer.parseInt(array[0]) > 0) {
            genotype = genotype + "A";
        }

        if (Integer.parseInt(array[1]) > 0) {
            genotype = genotype + "B";
        }

        if (genotype.length() == 0) {
            genotype = "N";
        }
        return genotype;
    }
}
