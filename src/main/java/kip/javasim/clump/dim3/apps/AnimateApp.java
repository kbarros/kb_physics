package kip.javasim.clump.dim3.apps;

import static scikit.util.Utilities.format;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.GrayScale;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Array3d;
import scikit.util.FileUtil;

public class AnimateApp extends Simulation {
	Grid3D grid = new Grid3D("Field");
	double time;
	Array3d a3d = new Array3d();
	
	public static void main(String[] args) {
		new Control(new AnimateApp(), "Clump Animation").getJob().throttleAnimation(true);
	}
	
	public void load(Control c) {
		c.frame(grid);
		params.add("Input directory", new DirectoryValue("/Users/kbarros/Desktop/Generated Data"));
		params.add("t start", 1.5);
		params.add("t finish", 700.0);
		params.add("dt", 20);
		params.add("lo", 0.0);
		params.add("hi", 10.);
		params.add("time");
	}

	public void animate() {
		params.set("time", format(time));
		grid.registerData(a3d);
	}
	
	public void clear() {
		grid.clear();
		params.set("time", 0);
	}
	
	public void run() {
		String indir = params.sget("Input directory");
		File outdir = FileUtil.getEmptyDirectory(indir, "images");
		DecimalFormat fmt = new DecimalFormat("0000");
		
		double lo = params.fget("lo");
		double hi = params.fget("hi");
		double ti = params.fget("t start");
		double tf = params.fget("t finish");
		double dt = params.fget("dt");
		
		grid.includeBoundary(false);
		grid.setColors(new GrayScale());
		grid.setScale(lo, hi);
//		grid.getSliceView().sliceResolution *= 2;
		
		for (time = ti; time < tf; time += dt) {
			String fname = indir+File.separator+"t="+format(time);
			a3d.readFile(new File(fname));
			
			Job.animate();
			try {
				String imgName = outdir+File.separator+fmt.format(time)+".png";
				ImageIO.write(grid.getImage(320, 320), "png", new File(imgName));
			} catch (IOException e) {
				System.err.println(e.getMessage());
			}
			Job.animate();
		}
		Job.animate();
	}
}
