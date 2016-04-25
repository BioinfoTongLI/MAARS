package org.micromanager.resultSaver;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.micromanager.cellstateanalysis.Cell;
import org.micromanager.cellstateanalysis.SpotsContainer;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.io.TmXmlWriter;

public class MAARSSpotsSaver {
	private String spotsXmlDir;
	private SpotsContainer container;

	public MAARSSpotsSaver(String pathToFluoDir) {
		spotsXmlDir = pathToFluoDir + "/spots/";
		if (!new File(spotsXmlDir).exists()) {
			new File(spotsXmlDir).mkdirs();
		}
	}

	public void saveSpots(String channel, SpotCollection spotsInChannel, String cellNb) {
		Model trackmateModel = container.getTrackmateModel();
		// for each cell
		File newFile = new File(spotsXmlDir + String.valueOf(cellNb) + "_" + channel + ".xml");
		TmXmlWriter spotsWriter = new TmXmlWriter(newFile);
		trackmateModel.setSpots(spotsInChannel, false);
		spotsWriter.appendModel(trackmateModel);
		try {
			spotsWriter.writeToFile();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void save(Cell cell) {
		this.container = cell.getSpotContainer();
		for (String channel : container.getUsingChannels()) {
			saveSpots(channel, container.getSpots(channel), String.valueOf(cell.getCellNumber()));
		}
	}
}