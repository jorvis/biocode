Introduction
============
We want this project to be made up of only well-documented scripts of general use. One-off 
scripts, or those with inappropriate hard-coded values (like databases, molecule values, 
etc.) should not be commited here. Previous collections I've seen became a mess because 
anything was allowed to go in it and no conventions were applied. The scripts weren't 
versioned, so it was impossible to keep track of how they changed over time. Also, contact 
information wasn't always included so it was difficult know who to suggest changes to or 
to inform if was broken.

We can make this Biocode repository a better place. By enforcing a few simple requirements 
and conventions it can shine like a beacon above many other trash-dump script repositories 
that have come before it. (drama.)

Details
=======

Before you add your script to this repository you'll need to do the following:

1. *Fork the repository* and learn how to send a pull request once you've written your 
   script(s).  You can learn more about that via [Github's documentation here](https://help.github.com/articles/using-pull-requests)

2. *Write documentation*.  Write it.  This should be embedded in your code in a way that 
   follows convention with its programming language. In Perl, for example, use perldoc. For 
   Python, use pydoc, etc. If you need guidance with this try starting with one of the template 
   scripts (found in the root directory) for your language as a reference (Perl and Python ones
   are provided.)  Suggested headings are:
   
   - NAME
   - SYNOPSIS
   - OPTIONS
   - DESCRIPTION
   - INPUT
   - OUTPUT
   - CONTACT

3. *Handle arguments well*.  Accept arguments with Getopt::Long (for perl), argparse (python) 
   or something equivalent for your language. Again, if you're not familiar with these you can 
   try using existing scripts as a template.

4. *Script naming*.  Name your scripts in all lower case to make everyone's lives easier. Too 
   often errors have been made because the user didn't know if the script was called FindOrfs.pl 
   or findORFs.pl or FiNDoRfS.pl or find_orfs.pl . There may not be one best convention, but it 
   is better to pick one and go with it. Use all lower case letters with an underscore between 
   words. Try to name your scripts with a long enough name to be descriptive, usually beginning 
   with a verb. For example, a name like 'get_seqs.pl' is crappy compared a longer, more 
   descriptive name like 'extract_assembly_sequences_from_chado.pl'.

5. *Respect return values*. The success or failure of a script is judged not by any log files it 
   might generate, but most directly by its return value. The Unix convention is that if a script 
   returns 0 it was successful, and any other value indicates failure. In perl, you can control 
   this by what you pass to the exit() function. calling exit() or exit(0) suggests the script 
   execution was peachy, while any other value, such as exit(42), indicates a failure. Any script 
   that fails but still returns 0 will wreak havoc on automated systems that may run it.
