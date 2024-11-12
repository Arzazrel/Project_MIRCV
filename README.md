# Project_MIRCV
**DESCRIPTION**  
  
Folder for the Multimedia Information Retrieval and Computer Vision project at the University of Pisa (year 2022-2023).  
This project was developed using IntelliJ IDEA 2023.2 (Community Edition) with Oracle OpenJDK version 19.0.2.  
This project consists of a local search engine developed completely in java that can operate on a collection of documents locally.  
Given a correctly formatted collection of documents, the program will create an inverted index, and other supporting data structures, which will allow textual queries, both conjunctive and disjunctive, to be performed on that collection.  
The program will implement several features to improve the information retrieval in the collection. These include: stopword removal, stemming, compression (of DIDs and term frequencies of posting lists), skipping with implementation of the MaxScore as a dynamic pruning algorithm.   
In addition, several testing features are implemented such as:  
- execution time tests can be performed on test queries passed in with correctly formatted collections;
- queries entered by the user can be repeated several times to test the average execution speed;
- various statistics may be calculated and displayed concerning the document collection (such as number of documents, average length of documents, number of poorly formatted documents, etc.);
- many of the techniques used can be tested individually (such as TUB calculation per term, compression and decompression testing and analysis, etc.);
- the size of all data structures saved on disk by the program can be shown;
- read tests can be carried out on the data structures saved on disk, which will show the data and read times.


**COLLECTION LINKS**  
  
Some collection, usefull for the operation of the search engine, can be found on the page accessible via this link [*msmarco/TREC-Deep-Learning-2020*](https://microsoft.github.io/msmarco/TREC-Deep-Learning-2020).
The collections used in the project are:
- test queries (2020) at this link [*test_2020_download*](https://msmarco.z22.web.core.windows.net/msmarcoranking/msmarco-test2020-queries.tsv.gz)
- test queries (2019) at this link [*test_2019_download*](https://msmarco.z22.web.core.windows.net/msmarcoranking/msmarco-test2019-queries.tsv.gz)
- document collection at this link [*doc_collection_download*](https://msmarco.z22.web.core.windows.net/msmarcoranking/collection.tar.gz)
  
The test queries collections must be placed in the '\Search_engine\src\main\resources\test' and will be used to test the performance of the search engine. 
These files are already present in the corresponding folder, so it should not be necessary to put them there.
The collection of documents must be placed in the '\Search_engine\src\main\resources' and will be used to create the inverted index and all the data structures required to do the
information retrieval in that document collection.
  
The names of the collections must not be changed. If you want to change them, you will have to change the name in the path in the ‘Costants’ class of the project, otherwise the 
collections will not be found by the programme and some functionality cannot be implemented.

**IDE Download**  
  
The version used in the project, like others, can be downloaded at this link [*Intellij_download*](https://www.jetbrains.com/idea/download/other.html). 


**JDK Download**  
  
To download Oracle JDK version 19, you can visit Oracle's official Java SE Downloads page at this link [*JDK19_download_page*](https://www.oracle.com/java/technologies/javase/jdk19-archive-downloads.html). 
Here, Oracle provides archived versions of JDK 19, inluded the 19.0.2 used for the project.
To access these downloads, you may need to accept the Oracle Technology Network License Agreement for Java SE, which typically allows for personal and development use but may 
have limitations for production environments without a commercial license. Once on the download page, look for "JDK 19.0.2" under the archived versions section, where you’ll find 
installers compatible with Linux, macOS, and Windows platforms.


**The folder contains:**  
  
- Search_engine: the folder containing the entire search engine java project
- Documentation: documentation regarding the application and the analysis carried out for the project.


**Project Import Guide**  
  
Note: you must not change the names of the files used, otherwise the prorgam cannot find them and cannot function correctly.  
  
0) Download the ‘Search_engine’ folder from the github repository.
1) Download the collection from the above link and put it in ‘Search_engine\src\main\resources’. 
   It is not necessary to unzip the collection, do not change the name of the collection otherwise the programme will not be able to find it.
2) If not present, download the test query collections from the links above and place them in ‘Search_engine\src\main\resources\test’.
3) Open IntelliJ IDEA.
4) Select 'File' from the top menu and then click on 'Open'.
5) Navigate to the folder containing the project ('Search_engine' folder), select it, and click OK (or Open).
6) IntelliJ will ask if you want to open the project as a new window or replace an existing one. Choose your preferred option.
5) Wait for IntelliJ to index and configure the project. If prompted to import the project configuration (especially for Maven or Gradle projects), choose the option to import.
6) Ensure that the correct JDK is selected:
     - Go to File -> Project Structure.
     - In the Project section, confirm the Project SDK version (program developed with version 19.0.2).
7) Build the Project by selecting 'Build' -> 'Build Project' from the menu.
8) To run the project, locate the main class or specific configuration you want to execute. Right-click on the main class file or Run Configuration and choose Run. The main class in this project is the class 'Main'.
   The project should now be running. If any dependencies are missing, IntelliJ might prompt you to download them automatically.


**Project Execution Guide**  
  
Example of basic operations.
0) Run the program.
1) A textual interface will be shown with various options divided into 3 categories:
      - Functionalities: includes the main functionalities of the search engine;
      - Statistics: includes the options to show the general statistics of the collection and data structures of the inverted index; 
      - Test: includes options for testing and verifying all the various search engine functionalities.
   To select an option, enter the corresponding character and press enter.
2) Press ‘i’ and send and then choose options for the various flags by typing ‘y’ to enable and ‘n’ to disable the corresponding flag. 
   Examples of settings:
   - Base case: 
      - is stopwords removal enabled : false
      - is compression enabled       : false
      - is scoring BM25 enabled      : false
      - is skipping enabled          : false
      - is dynamic pruning enabled   : false
      - delete partial file enabled  : false
      - whole PL in memory enabled   : false
    - Only Stopwords removal:
      - is stopwords removal enabled : true
      - is compression enabled       : false
      - is scoring BM25 enabled      : false
      - is skipping enabled          : false
      - is dynamic pruning enabled   : false
      - delete partial file enabled  : false
      - whole PL in memory enabled   : false 
    - All functionality enabled:
      - is stopwords removal enabled : true
      - is compression enabled       : true
      - is scoring BM25 enabled      : false
      - is skipping enabled          : true
      - is dynamic pruning enabled   : true
      - delete partial file enabled  : true
      - whole PL in memory enabled   : true
  
    Note: 
      - Partial inverted indexes, contained in the ‘partial’ folder, are no longer used once the complete inverted index is finished, so it is recommended to save space to enable 
        the flag to delete them once the inverted index is finished. The presence or absence of these partial files will not affect the operation of the search engine, 
	it will only change the amount of disk used by the programme.
      - The scoring flag only changes the type of scoring function used: if ‘true’ it uses BM25 and if ‘false’ it uses TFIDF.
3) Wait for the program to finish creating the inverted index and all structures necessary for operation. This can take more than 30 minutes.
4) Once finished, the text interface will be displayed again.
5) Try out as many queries as you like.  
   For each query, follow this procedure:
   - Type ‘q’ and enter.
   - Enter the query you want and press enter.
   - Choose the query type: ‘d’ for disjunctive and ‘c’ for conjunctive and press enter.
   - Choose the number of results you want to see: 10 or 20.
   - The ordered results and the query execution time will be displayed.
6) Close the programme by typing ‘x’ and pressing enter.

Tips:  
   - Please note that you can also change flags with f but depending on the change it may take several minutes and even the whole inverted index rebuild. 
     It is advisable to use this functionality only for the flags ‘dynamic pruning’, ‘delete partial file’ and ‘whole PL in memory’, which make non-burdensome changes. 
     Even changing the scoring function can be done without too much trouble, due to the recalculation of the Term Upper Bounds, which should not take more than a few minutes.
     For all other flags changes, it is recommended to exit the program with ‘x’, then re-run the program and select ‘i’ and enter the new flag values. 
   - If once the inverted index has been made, it is desired to make another one with new flag values, it is recommended to exit the program with ‘x’, then re-run the program and select ‘i’
     and enter the new flag values.
   - If you want to keep several flags settings so that you can work with them without having to rebuild the inverted index each time, I recommend saving the contents of the ‘merged’ and 
     ‘upperBound’ folders and the ‘collectionStatistics’ and ‘flags’ files. Save these files of a configuration, to make the search engine run with those settings, simply copy those 
     folders and files into the ‘resources’ folder and then start the program. 
	    

**Developer's notes**  
  
The work related to the university examination has been done and the project is completed. 
There may be updates or improvements to the project in the future, but nothing is planned for now.


**Credits**  
  
This project was realised in collaboration with Martina Marino and Roberta Matrella, who were fundamental to the implementation and correct functioning of the search engine code when this repository was started. 


**Developers:**  
- Alessandro Diana