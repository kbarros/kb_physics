package kip.javasim.clump.dim2.apps;

import static scikit.util.Utilities.format;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.GrayScale;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.FileUtil;

public class AnimateApp extends Simulation {
	Grid grid = new Grid("Field");
	double time;

	public static void main(String[] args) {
		new Control(new AnimateApp(), "Clump Animation").getJob().throttleAnimation(true);
	}
	
	public void load(Control c) {
		c.frame(grid);
		params.add("Input directory", new DirectoryValue());
		params.add("t start", 1700.0);
		params.add("t finish", 2150.0);
		params.add("dt", 1);
		params.add("lo", 0.2);
		params.add("hi", 4.);
		params.add("time");
	}

	public void animate() {
		params.set("time", format(time));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		String indir = params.sget("Input directory");
		File outdir = FileUtil.getEmptyDirectory(indir, "images");
		DecimalFormat fmt = new DecimalFormat("0000");
		
		double[] data = null;
		double lo = params.fget("lo");
		double hi = params.fget("hi");
		double ti = params.fget("t start");
		double tf = params.fget("t finish");
		double dt = params.fget("dt");
		
		for (time = ti; time < tf; time += dt) {
			String fname = indir+File.separator+"t="+format(time);
			try {
				DataInputStream dis = FileUtil.disFromString(fname);
				int w = dis.readInt();
				int h = dis.readInt();
				if (data == null)
					data = new double[w*h];
				for (int i = 0; i < w*h; i++)
					data[i] = dis.readFloat();
				dis.close();
				
				grid.setColors(new GrayScale());
				grid.setScale(lo, hi);
				grid.registerData(w, h, data);
				try {
					String imgName = outdir+File.separator+fmt.format(time)+".png";
					ImageIO.write(grid.getImage(), "png", new File(imgName));
				} catch (IOException e) {
				}
			}
			catch (IOException e) {
				System.err.println("Can't read file: " + fname);
			}
			Job.animate();
		}
		Job.animate();
	}
}
