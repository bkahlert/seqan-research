#!/usr/bin/env python
"""Helper code for app tests.

This module contains helper functions and classes for making app tests easy.
The advantage of using Python for this is easier portability instead of relying
on Unix tools such as bash and diff which are harder to install on Windows than
Python.

App tests are performed by executing the programs on test data and comparing
their output to previously generated "golden" output files.

Classes/Functions:

  class TestConf -- stores configuration of a test.
  class TestPathHelper -- helps with constructing paths.
  function runTest -- runs a test configured by a TestConf object.
  function autolocateBinary -- locates a binary, possibly in an intermediary
                               directory.
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import difflib
import logging
import optparse
import os
import os.path
import re
import subprocess
import shutil
import sys
import tempfile


# Valgrind flags, taken from CMake output, ideally given to test script by CMake?
SUPPRESSIONS = '--suppressions=' + os.path.join(os.path.dirname(__file__), '..', '..', 'seqan.supp')
VALGRIND_FLAGS = [SUPPRESSIONS] + '--error-exitcode=1 -q --tool=memcheck --leak-check=yes --show-reachable=yes --workaround-gcc296-bugs=yes --num-callers=50 --'.split()
VALGRIND_PATH = '/group/agabi/software/bin/valgrind'

class TestConf(object):
    """Configuration for one tests.

    A test configuration consists of the parameters to give to the
    program and the expected result.

    Attrs:
      program -- string, path to binary to execute.
      args -- list of strings with arguments to the program.
      to_diff -- optional list of pairs with (output-file, expected-file) paths
                 diff, the contents of output-file should be equal to the
                 contents of expected-file.
      name -- optional string, name of the test.
      redir_stdout -- optional string that gives the path to redirect stdout to
                      if the variable is not None.
    """

    def __init__(self, program, args, to_diff=[], name=None,
                 redir_stdout=None):
        """Constructor, args correspond to attrs."""
        self.program = program
        self.args = args
        self.to_diff = to_diff
        self.name = name
        self.redir_stdout = redir_stdout
        if not hasattr(TestConf, 'valgrind'):
            self.valgrind = False
        else:
            self.valgrind = TestConf.valgrind

    def __str__(self):
        fmt = 'TestConf(%s, %s, %s, %s, %s)'
        return fmt % (repr(self.program), self.args, self.to_diff, self.name,
                      self.redir_stdout)


class TestPathHelper(object):
    """Helper class for paths.

    TestPathHelper objects are configured with the appropriate paths.  The
    provide functions to construct when executing tests.
    """

    def __init__(self, source_base_path, binary_base_path,
                 tests_dir):
        self.temp_dir = None
        self.source_base_path = source_base_path
        self.binary_base_path = binary_base_path
        self.tests_dir = tests_dir
        self.created_paths = []

    def inFile(self, path):
        """Convert the path of a test file.

        The given path, relative to the test directory, will be transformed into an
        absolute path to the file.

        Args:
          path -- relative path to a file in the test directory.

        Returns:
          Absolute to the file.
        """
        result = os.path.join(self.source_base_path, self.tests_dir, path)
        logging.debug('inFile(%s) = %s', path, result)
        return result

    def outFile(self, path, subdir=None):
        """Convert the path of an output file.

        The given path will be converted to a path to a temporary file.  The path
        to this file will be created.

        If subdir is set then a subdirectory with this name will be created and
        the output will be relative to subdir.
        """
        if not self.temp_dir:
            self.temp_dir = tempfile.mkdtemp()
            if not os.path.isdir(self.temp_dir):
                self.created_paths.append(self.temp_dir)
                os.makedirs(self.temp_dir)
        target_dir = self.temp_dir
        if subdir:
            target_dir = os.path.join(self.temp_dir, subdir)
            if not os.path.isdir(target_dir):
                self.created_paths.append(target_dir)
                os.makedirs(target_dir)
        logging.debug('outFile(%s, %s) = %s', path, subdir, self.temp_dir)
        res = os.path.join(target_dir, path)
        self.created_paths.append(res)
        return res

    def deleteTempDir(self):
        """Remove the temporary directory created earlier and all files below."""
        print >>sys.stderr, 'DELETING TEMP DIR', self.temp_dir
        if self.temp_dir:
            shutil.rmtree(self.temp_dir)


def autolocateBinary(base_path, relative_path, binary_name):
  """Autolocates a binary, possibly in an intermediary path.

  When building applications with CMake, they do not always have the same
  relative path from the binary build directory.  For Unix Makefiles, the path
  could be 'core/apps/tree_recon' whereas for Visual Studio, it could be
  'core/apps/Release/tree_recon'.

  This function tries to automatically guess the name of the file and return
  the first one it finds.
  """
  # Names of intermediary directories and possible file extensions.
  intermediary_dir_names = ['', 'Debug', 'Release']
  extensions = ['', '.exe']
  paths = []
  # Try all possible paths.
  for dir_name in intermediary_dir_names:
    for ext in extensions:
      res_list = [base_path, relative_path, dir_name, binary_name + ext]
      filtered_list = [x for x in res_list if x]  # Filter out empty strings.
      paths.append(os.path.join(*filtered_list))
      if dir_name:
        paths.append('/'.join([base_path] + relative_path.split('/')[:-1] +
                              [dir_name] + relative_path.split('/')[-1:] +
                              [binary_name]))
  for path in paths:
    logging.debug('Trying path %s', path)
    if os.path.isfile(path):
      logging.debug('  Found binary %s', path)
      return path
  # Fall back ot Unix default.
  return os.path.join(base_path, relative_path, binary_name)


def runTest(test_conf):
    """Run the test configured in test_conf.

    Args:
      test_conf -- TestConf object to run test for.

    Returns:
      True on success, False on any errors.

    Side Effects:
      Errors are printed to stderr.
    """
    # Execute the program.
    logging.debug('runTest(%s)', test_conf)
    args = [test_conf.program] + test_conf.args
    if test_conf.valgrind:
        # Call through valgrind.
        args = [VALGRIND_PATH] + VALGRIND_FLAGS + args
    logging.debug('Executing "%s"', ' '.join(args))
    stdout_file = subprocess.PIPE
    if test_conf.redir_stdout:
        logging.debug('  Redirecting stdout to "%s".' % test_conf.redir_stdout)
        stdout_file = open(test_conf.redir_stdout, 'w+')
    try:
        process = subprocess.Popen(args, stdout=stdout_file,
                                   stderr=subprocess.PIPE)
        retcode = process.wait()
        logging.debug('  return code is %d', retcode)
        if retcode != 0:
            fmt = 'Return code of command "%s" was %d.'
            print >>sys.stderr, '--- stdout begin --'
            print >>sys.stderr, fmt % (' '.join(args), retcode)
            print >>sys.stderr, stdout_file.read()
            print >>sys.stderr, '--- stdout end --'
            stdout_file.close()
            stderr_contents = process.stderr.read()
            print >>sys.stderr, '-- stderr begin --'
            print >>sys.stderr, stderr_contents
            print >>sys.stderr, '-- stderr end --'
            return False
    except Exception, e:
        fmt = 'ERROR (when executing "%s"): %s'
        if stdout_file is not subprocess.PIPE:
            stdout_file.close()
        print >>sys.stderr, fmt % (' '.join(args), e)
        return False
    # Handle error of program, indicated by return code != 0.
    if retcode != 0:
        print >>sys.stderr, 'Error when executing "%s".' % ' '.join(args)
        print >>sys.stderr, 'Return code is %d' % retcode
        if stdout_file is not subprocess.PIPE:
            stdout_file.seek(0)
        stdout_contents = process.stdout.read()
        if stdout_contents:
            print >>sys.stderr, '-- stdout begin --'
            print >>sys.stderr, stdout_contents
            print >>sys.stderr, '-- stdout end --'
        else:
            print >>sys.stderr, '-- stdout is empty --'
        stderr_contents = process.stderr.read()
        if stderr_contents:
            print >>sys.stderr, '-- stderr begin --'
            print >>sys.stderr, stderr_contents
            print >>sys.stderr, '-- stderr end --'
        else:
            print >>sys.stderr, '-- stderr is empty --'
    # Close standard out file if necessary.
    if stdout_file is not subprocess.PIPE:
        stdout_file.close()
    # Compare results with expected results, if the expected and actual result
    # are not equal then print diffs.
    result = True
    for tuple_ in test_conf.to_diff:
        expected_path, result_path = tuple_[:2]
        transforms = []
        if len(tuple_) >= 3:
            transforms += tuple_[2]
        try:
            with open(expected_path, 'r') as f:
                expected_str = f.read()
            for t in transforms:
                expected_str = t.apply(expected_str, True)
            with open(result_path, 'r') as f:
                result_str = f.read()
            for t in transforms:
                result_str = t.apply(result_str, False)
            if expected_str == result_str:
                continue
            fmt = 'Comparing %s against %s'
            print >>sys.stderr, fmt % (expected_path, result_path)
            diff = difflib.unified_diff(expected_str.splitlines(),
                                        result_str.splitlines())
            for line in diff:
                print >>sys.stderr, line
            result = False
        except Exception, e:
            fmt = 'Error when trying to compare %s to %s: %s ' + str(type(e))
            print >>sys.stderr, fmt % (expected_path, result_path, e)
            result = False
    return result


class ReplaceTransform(object):
    """Transformation on left and/or right files to diff."""

    def __init__(self, needle, replacement, left=False, right=False):
        self.needle = needle
        self.replacement = replacement
        self.left = left
        self.right = right

    def apply(self, text, is_left):
        if (is_left and not self.left) or (not is_left and not self.right):
            return text  # Skip if no transform is to be applied.
        return text.replace(self.needle, self.replacement)


class NormalizeScientificExponentsTransform(object):
	"""Transformation that normalized scientific notation exponents.
	
	On Windows, scientific numbers are printed with an exponent padded to
	a width of three with zeros, e.g. 1e003 instead of 1e03 as on Unix.
	
	This transform normalizes to Unix or Windows.
	"""
	
	def __init__(self, normalize_to_unix=True):
		self.normalize_to_unix = normalize_to_unix
	
	def apply(self, text, is_left):
		"""Apply the transform."""
		if self.normalize_to_unix:
			return re.sub(r'([-+]?[0-9]*\.[0-9]+[eE][-+]?)0([0-9]{2})', r'\1\2', text)
		else:
			return re.sub(r'([-+]?[0-9]*\.[0-9]+[eE][-+]?)([0-9]{2})', r'\10\2', text)


class RegexpReplaceTransform(object):
    """Transformation that applies regular expression replacement."""
    
    def __init__(self, needle, replacement, left=True, right=True):
        self.needle = needle
        self.replacement = replacement
        self.left = left
        self.right = right
    
    def apply(self, text, is_left):
        """Apply the transform."""
        if (is_left and not self.left) or (not is_left and not self.right):
            return text  # Skip if no transform is to be applied.
            return re.sub(self.needle, self.replacement, text)


def main(main_func):
    """Run main_func with the first and second positional parameter.""" 
    parser = optparse.OptionParser("usage: run_tests [options] SOURCE_ROOT_PATH BINARY_ROOT_PATH")
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true')
    parser.add_option('--valgrind', dest='valgrind', action='store_true')
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error('Incorrect number of arguments!')
        return 2
    if options.verbose:
        logging.root.setLevel(logging.DEBUG)
    if options.valgrind:
        TestConf.valgrind = True
    return main_func(args[0], args[1])
