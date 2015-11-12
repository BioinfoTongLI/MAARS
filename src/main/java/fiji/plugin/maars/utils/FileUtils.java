package fiji.plugin.maars.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.micromanager.internal.utils.ReportingUtils;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.io.TmXmlWriter;
import ij.IJ;

/**
 * @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 4, 2015
 */
public class FileUtils {

	public static boolean isValid(String path) {
		File f = new File(path);
		if (f.exists()) {
			return true;
		} else {
			return false;
		}
	}
	
	public static ArrayList<File> listFilesForFolder( File folder) {
		ArrayList<File> list = new ArrayList<File>();
	    for ( File fileEntry : folder.listFiles()) {
	        if (fileEntry.isDirectory()) {
	            listFilesForFolder(fileEntry);
	        } else {
	            list.add(fileEntry);
	        }
	    }
	    return list;
	}
	
	public static void writeSpotFeatures(File file, Model model){
		File newFile = new File(file.getParentFile() + "/" + file.getName() +".xml");
		ReportingUtils.logMessage("Writing to :" + newFile.toString());
		TmXmlWriter writer = new TmXmlWriter( newFile , model.getLogger() );
		writer.appendModel( model );
		try {
			writer.writeToFile();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Convert an unix path in windows path if program is running on windows OS
	 * @param unixPath
	 * @return String path
	 */
	public static String convertPath(String unixPath) {
		String path = unixPath;
		ReportingUtils.logMessage("program running on windows : "+IJ.isWindows());
		ReportingUtils.logMessage("path is containing '/' : "+path.contains("/"));
		if (IJ.isWindows() && path.contains("/")) {
			path = path.replace("/", "\\\\");
		}
		return path;
	}
}
