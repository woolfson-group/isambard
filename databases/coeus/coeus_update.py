import datetime
import logging
import os
import pickle
import shutil

from add_ons.filesystem import pdbe_status_code, current_codes_from_pdb
from settings import global_settings


# TODO Re-write this file for fitting in with coeus structure.
# TODO Need to deal with not having Bio.PDB.PDBList.


def get_retry_codes(database_folder=global_settings["interactions_database"]["folder"]):
    """ Get list of codes that did not have PDBE files at the previous attempt

    Many structures that are new to the PDB are not yet present in the PDBE.
    The codes for such structures are returned by this function so that they can be retried on the next update.

    Parameters
    ----------
    database_folder : str, optional
        Parent folder containing PDB-style structure and PDB files.
        Defaults to database folder in isambard_global_settings.

    Returns
    -------
    try_agains : list
        List of all PDB codes in PDB that were not yet in PDBE at last update attempt.
    """
    codes = current_codes_from_pdb()
    try_agains = []
    for code in codes:
        p = "{0}/{1}/{2}/structures/".format(database_folder, code[1:3], code)
        if len(os.listdir(p)) == 0:
            try_agains.append(code)
    return try_agains


class DBUpdate:
    """ For updating the local PDB directory, as well as the database directory and its files.

    Uses Bio.PDB.PDBList to get the three lists (new, modified and obsolete) of PDB codes from the latest version
    the PDB on the web. These lists then are each attributes of the class.
    A log file is created, and stored in the 'update_logs' folder within the database directory.
    The log file has time stamps for each set of changes made.
    The first update of the day creates a new log file named after that date (yyyy-mm-dd_update.log).
    If multiple updates are run on a single day (unnecessary), then these logs will be appended to the same file.
    The functions of the class update the PDB directory, the database directory, and the sql database.

    Attributes
    ----------
    database : str
        Database folder in isambard_global_settings.
    today : datetime.date
        Today's date, used for logging.
    logfolder : str
        Folder for storing log file.
        Stored within database folder, subfolder named according to the year of the update.
    logfile : str
        Filepath for the logfile for the update.
    logger : logger.Logger
        Logging object for constructing the logfile and its contents.
    fh : logger.FileHandler
        File handle for the logfile.
    obsolete_folder : str
        Subfolder in the database folder for storing obsolete PDB files.
    pl : Bio.PDB.PDBList.PDBList
    list_of_pdbs : list of list of str
        Compiled new, modified and obsolete PDB codes from the PDB.
    new : list of str
        Recently added PDB codes in the PDB.
    modified : list of str
        Recently modified PDB codes in the PDB.
    obsolete : list of str
        PDB codes in the PDB that have recently become obsolete.
    """
    def __init__(self):

        # Set up everything required for the logging.
        self.database = global_settings["structural_database"]["path"]
        self.today = datetime.date.today()
        self.logfolder = "{0}update_logs/{1}".format(self.database, self.today.year)
        if not os.path.exists(self.logfolder):
            os.mkdir(self.logfolder)
        # Create logfile named after the date on which DBUpdater is run.
        self.logfile = ("{0}/{1}_update.log".format(self.logfolder, self.today.isoformat()))
        self.logger = logging.getLogger("{0}.DBUpdate".format(__name__))
        self.logger.setLevel(logging.DEBUG)
        self.fh = logging.FileHandler(self.logfile)
        self.fh.setLevel(logging.DEBUG)
        self.logger.addHandler(self.fh)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.fh.setFormatter(formatter)

        self.obsolete_folder = "{0}obsolete".format(self.database)
        if not os.path.exists(self.obsolete_folder):
            os.mkdir(self.obsolete_folder)

        self.pl = PDBList()
        self.list_of_pdbs = self.pl.get_recent_changes()

        (new, modified, obsolete) = self.list_of_pdbs
        # Safe-keeping: check all entries are indeed PDB codes.
        self.new = [code.lower() for code in new if len(code) == 4]
        self.modified = [code.lower() for code in modified if len(code) == 4]
        self.obsolete = [code.lower() for code in obsolete if len(code) == 4]

    def pickle_pdb_lists(self):
        """ Pickle the new, modified and obsolete PDB codes as a dictionary.

        Compiles a dictionary and stores it as a pickle file within the logfolder.
        Keys are 'new', 'modified' and 'obsolete'.
        Values are lists of the associated PDB codes.

        Returns
        -------
        None
        """
        pickle_file = "{0}/{1}_pdb_lists.p".format(self.logfolder, self.today.isoformat())
        p_dict = {'new': self.new,
                  'modified': self.modified,
                  'obsolete': self.obsolete
                  }
        pickle.dump(p_dict, open(pickle_file, 'wb'))
        return

    def run_update(self):
        """ Update all database files for each PDB code from the PDB list.

        In each case, the 'files' are all files associated with the FileSystem class object.
        Create files for new PDB codes.
        Delete and re-create files for modified PDB codes.
        Move files for obsolete PDB codes into the folder database_folder/obsolete.

        Returns
        -------
        None
        """
        retry_file = '{0}retry_codes/retry_code_list.p'.format(self.database)
        if not os.path.exists(retry_file):
            raise IOError('{0} does not exist. Should be a pickle file of retry codes'.format(retry_file))
        retries = pickle.load(open(retry_file, 'rb'))
        for code in retries:
            if pdbe_status_code(code) != 404:
                self.new.append(code)
                retries.remove(code)

        # Create files for each code in new.
        for code in self.new:
            if pdbe_status_code(code) == 404:
                retries.append(code)
                self.logger.info("Added to retries list new: {0}".format(code))
                continue
            try:
                dbi = DatabaseItem(code=code)
                fill = populate_coeus.fill_database(dbi=dbi)
                if fill:
                    self.logger.info("Successfully populated prot_graph with new {0}".format(code))
                else:
                    self.logger.debug("Could not populate prot_graph with new {0}".format(code))
                    continue
            except ValueError as e:
                self.logger.debug("Problem in socket file for new {0}\n{1}".format(code, e))
            except:
                self.logger.debug("Could not get DatabaseItem for new {0}".format(code))
                continue

        for code in self.modified:
            # Remove data from the database for code
            Pdb.objects.filter(pdb=code).delete()
            self.logger.info("Removed old data from prot_graph for modified {0}".format(code))

            # Remove existing files
            code_folder = "{0}{1}/{2}".format(self.database, code[1:3], code)
            if os.path.exists(code_folder):
                shutil.rmtree(code_folder)
            self.logger.info("Removed old files from db folder for modified {0}".format(code))

            # Replace with new files
            if pdbe_status_code(code) == 404:
                retries.append(code)
                self.logger.info("Added to retries list modified: {0}".format(code))
                continue
            try:
                dbi = DatabaseItem(code=code)

                if len(dbi.filesystem.mmols) == 0:
                    self.logger.debug("Could not download mmol files for modified {0}".format(code))

                fill = populate_coeus.fill_database(dbi=dbi)
                if fill:
                    self.logger.info("Successfully populated prot_graph with modified {0}".format(code))
                else:
                    self.logger.debug("Could not populate prot_graph with modified {0}".format(code))
                continue

            except:
                self.logger.debug("Could not get DatabaseItem for modified {0}".format(code))
                continue

        for code in self.obsolete:
            # Remove data from database for code
            Pdb.objects.filter(pdb=code).delete()
            self.logger.info("Removed old data from prot_graph for obsolete {0}".format(code))

            # Move files for each code in obsolete to the 'obsolete' subdirectory in database_dir.
            code_folder = "{0}{1}/{2}".format(self.database, code[1:3], code)
            if os.path.exists(code_folder):
                destination_code_folder = "{0}/{1}/{2}".format(self.obsolete_folder, code[1:3], code)

                if os.path.exists(destination_code_folder):
                    shutil.rmtree(destination_code_folder)

                shutil.move(code_folder, destination_code_folder)
                self.logger.info("Moved {0} to {1}".format(code_folder, destination_code_folder))

                # Clean up if there are no other codes with the same middle two letters.
                two_letter_folder = "{0}{1}".format(self.database, code[1:3])
                if len(os.listdir(two_letter_folder)) == 0 and two_letter_folder != self.database:
                    os.rmdir(two_letter_folder)

        retries = list(set(retries))
        pickle.dump(retries, open(retry_file, 'wb'))

        return


__author__ = 'Jack W. Heal'
__status__ = 'Development'
