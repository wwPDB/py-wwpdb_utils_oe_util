##
# File:    OeDescriptorUtils.py
# Author:  jdw
# Date:    22-May-2017
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to standardize chemical descriptors

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import logging

from openeye import oechem

logger = logging.getLogger(__name__)


class OeDescriptorUtils(object):
    ''' Utilities to standardize chemical descriptors
    '''

    def __init__(self):
        pass
        #

    def standardizeSmiles(self, smiles, type="ISOMERIC"):  # pylint: disable=redefined-builtin
        """ Return a standardized SMILES (type) or None
        """
        smilesOut = None
        try:
            mol = oechem.OEGraphMol()
            if (oechem.OEParseSmiles(mol, smiles) == 1):
                oechem.OEAssignAromaticFlags(mol)
                if type == "CANNONICAL":
                    smilesOut = oechem.OECreateCanSmiString(mol)
                elif type == "ISOMERIC":
                    smilesOut = oechem.OECreateIsoSmiString(mol)
            else:
                logger.error("Unable to parse input SMILES '%s'", smiles)

        except Exception as e:
            logger.exception("Error '%s' occured. Arguments %s.", str(e), e.args)

        return smilesOut
