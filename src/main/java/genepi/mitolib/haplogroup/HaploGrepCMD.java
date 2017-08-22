package genepi.mitolib.haplogroup;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import genepi.base.Tool;
import server.main.HaplogrepCMD;


public class HaploGrepCMD  extends Tool {

	String[] args;
	
	public String[] getArgs() {
		return args;
	}


	public void setArgs(String[] args) {
		this.args = args;
	}


	public HaploGrepCMD(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("perform a haplogroup classification based on HaploGrep 2\n\n");
	}

	@Override
	public void createParameters() {

	}

	
	
	@Override
	public int run() {
	
		HaplogrepCMD.main(args); 
	
		return 0;
	}
	
	
	}


