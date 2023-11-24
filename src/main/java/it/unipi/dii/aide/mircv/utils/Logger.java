package it.unipi.dii.aide.mircv.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

import static it.unipi.dii.aide.mircv.utils.Constants.DEBUG_FOLDER;

public final class Logger {

    private static String logFileName;

    public static Logger dict_logger = new Logger("dict.txt");
    public static Logger docId_logger = new Logger("docid.txt");
    public static Logger termFreq_logger = new Logger("termFreq.txt");
    public static Logger collStats_logger = new Logger("collStats.txt");
    public static Logger docTable_logger = new Logger("docTable.txt");

    public Logger(String logFileName) {
        Logger.logFileName = DEBUG_FOLDER + logFileName;
    }

    public void logInfo(String message) {
        log("INFO", message);
    }

    public void logWarning(String message) {
        log("WARNING", message);
    }

    public void logError(String message) {
        log("ERROR", message);
    }

    private void log(String logLevel, String message) {
        String formattedLog = getFormattedLog(logLevel, message);

        try (PrintWriter writer = new PrintWriter(new FileWriter(logFileName, true))) {
            writer.println(formattedLog);
        } catch (IOException e) {
            System.err.println("Failed to write to log file: " + e.getMessage());
        }
    }

    private String getFormattedLog(String logLevel, String message) {
        Date now = new Date();
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String formattedDate = dateFormat.format(now);
        return "[" + formattedDate + "] [" + logLevel + "] " + message;
    }
}
