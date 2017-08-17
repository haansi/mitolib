package genepi.ngs;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;

import org.apache.commons.compress.utils.IOUtils;
import org.junit.Test;

import genepi.bam.VariantBuilder;
import junit.framework.Assert;

public class TestBam2Var {

	@Test
	public void testHG01500_low_coverage() {
		
		String in = "data/HG01500.IBS.low_coverage.MT.bam";
		String out = "output";
		String ref = "data/rcrs.fasta";
		double vaf = 0.1;
		int qual = 30;
		int mapqual = 200;

		VariantBuilder builder = new VariantBuilder(in);
		builder.setReference(ref);
		builder.setOutDirectory(out);
		builder.setVaf(vaf);
		builder.setQual(qual);
		builder.setMapqual(mapqual);
		try {
			builder.build();
		} catch (MalformedURLException e) {
			fail("URL failed");
			e.printStackTrace();
	
		} catch (IOException e) {
			fail("IO failed");
			e.printStackTrace();
		
		}
		String test="";

        try {
			test = new java.util.Scanner(new File("output/HG01500.IBS.low_coverage.MT.bam.txt.hsd"),"UTF8").useDelimiter("\\n").next();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        System.out.println(test);
       assertTrue(test.equals("HG01500.IBS.low_coverage.MT.bam	1-16569		146C	263G	750G	1438G	4769G	6776C	8860G	12957C	15326G	16519C	"));
	}

	
	@Test
	public void testHG01500_exome() {
		
		String in = "data/HG01500.IBS.low_coverage.MT.bam";
		String out = "output";
		String ref = "data/rcrs.fasta";
		double vaf = 0.1;
		int qual = 30;
		int mapqual = 200;

		VariantBuilder builder = new VariantBuilder(in);
		builder.setReference(ref);
		builder.setOutDirectory(out);
		builder.setVaf(vaf);
		builder.setQual(qual);
		builder.setMapqual(mapqual);
		try {
			builder.build();
		} catch (MalformedURLException e) {
			fail("URL failed");
			e.printStackTrace();
	
		} catch (IOException e) {
			fail("IO failed");
			e.printStackTrace();
		
		}
		String test="";

        try {
			test = new java.util.Scanner(new File("output/HG01500.IBS.exome.MT.bam.txt.hsd"),"UTF8").useDelimiter("\\n").next();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
       assertTrue(test.equals("HG01500.IBS.exome.MT.bam	1-16569		146C	263G	750G	1438G	3107T	4769G	6776C	8860G	12957C	15326G	16519C	"));
	}
	
	
}
