import os
import shutil
import subprocess
import struct

from biased_ff import IO  #CHANGE


class FragmentQualityAnalyzer:
    '''Analyze the fragment quality of a design'''
    def __init__(self, runpsipred_single, csblast_path,
            blastpgp, placeholder_seqs, sparksx_path, 
            fragment_picker_cmd, vall, weights_file, 
            rosetta_database=None):
        
        self.runpsipred_single_abs = os.path.abspath(runpsipred_single)
        self.csblast_path_abs = os.path.abspath(csblast_path)
        self.blastpgp_abs = os.path.abspath(blastpgp)
        self.placeholder_seqs_abs = os.path.abspath(placeholder_seqs)
        self.sparksx_path_abs = os.path.abspath(sparksx_path)

        self.fragment_picker_cmd = fragment_picker_cmd
        self.vall_abs = os.path.abspath(vall)
        self.weights_file_abs = os.path.abspath(weights_file)
        self.rosetta_database = rosetta_database
     
    def run_psipred(self, fasta_file):
        '''Run psipred to predict the secondary structures of a sequence.
        Return the name of the ss2 file.
        '''
        subprocess.check_call([self.runpsipred_single_abs, fasta_file])
        return os.path.basename(fasta_file)[:-5] + 'ss2'
     
    def run_csblast(self, fasta_file):
        '''Use CS-Blast to generate the sequence profile of a sequence.
        Return the name of the binary checkpoint file.
        '''
        checkpoint_binary = os.path.basename(fasta_file)[:-5] + 'check'
        
        subprocess.check_call([os.path.join(self.csblast_path_abs, 'bin', 'csbuild'),
                               '-i', fasta_file, '-I', 'fas', 
                               '-D', os.path.join(self.csblast_path_abs, 'data', 'K4000.crf'),
                               '-o', checkpoint_binary, '-O', 'chk'])

        return checkpoint_binary
   
    def parse_checkpoint_file(self, binary_file, text_file, sequence_length):
        '''Convert a binary checkpoint file to a text checkpoint file.'''
        column_map = ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18 )
        
        with open(binary_file, 'rb') as bf:
            data = bf.read()
            offset = 4 + sequence_length
            sequence = struct.unpack(str(sequence_length) + 'c', data[4:offset]) 

            with open(text_file, 'w') as tf:
                tf.write(str(sequence_length) + '\n')

                for i in range(sequence_length):
                    tf.write(sequence[i].decode())
                    profile = struct.unpack('20d', data[offset + i * 160 : offset + (i + 1) * 160])
                    
                    for j in column_map:
                        tf.write(' {:6.4f}'.format(profile[j]))
                    tf.write('\n')


                tf.write('END')
    
    def create_pssm(self, fasta_file):
        '''Create the pssm file.
        Return the name of the pssm file.
        '''
        blast_file = os.path.basename(fasta_file)[:-5] + 'blast'
        pssm_file = os.path.basename(fasta_file)[:-5] + 'pssm'
        sequence = IO.fasta_file_to_sequences(fasta_file)[1][0]


        with open(blast_file, 'w') as f:
            f.write(os.path.basename(fasta_file)[:-5]  + ' ')
            f.write(sequence)

        subprocess.check_call([self.blastpgp_abs, '-i', fasta_file, '-B', blast_file, '-Q', pssm_file, 
                               '-t', '1', '-j', '1', '-h', '0.001', '-e', '0.001', '-b', '0', '-k', '0',
                               '-d', self.placeholder_seqs_abs])

        return pssm_file

    def run_sparksx(self, fasta_file, pssm_file):
        '''Run the sparks-x program.
        Return the name of the phipsi file. 
        '''
        shutil.copy(pssm_file, fasta_file + '.pssm')
        exitcode, err = 0, None
        try:
            subprocess.check_call([os.path.join(self.sparksx_path_abs, 'bin', 'psiblast.sh'), fasta_file], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            exitcode, err = e.returncode, e.output
            print("x{} e{}".format(exitcode,err))
        subprocess.check_call([os.path.join(self.sparksx_path_abs, 'SPINE-X', 'spineX.pl'), fasta_file + '.pssm'])
        subprocess.check_call(['python2', os.path.join(self.sparksx_path_abs, 'bin', 'buildinp.py'), 
          fasta_file + '.phipsi', fasta_file + '.pssm', fasta_file + '.inp'])

        return fasta_file + '.phipsi'

    def clean(self, fasta_file):
        '''Clean all the temporary files that are no longer needed after fragment picking.'''
        base = os.path.basename(fasta_file)[:-5]
        files_to_clean = [base + ex for ex in
                ['blast', 'check', 'checkpoint', 'fasta.phipsi', 'fasta.pssm', 'fasta.inp',
                 'horiz', 'pssm', 'ss', 'ss2']] 
        
        for f in files_to_clean:
            if os.path.exists(f):
                os.remove(f)

    def pick_fragments(self, input_pdb, input_fasta, working_dir):
        '''Pick fragments of a structure.
        Return the name of the fragment describing file.
        '''
        cwd = os.getcwd()
        input_pdb_abs = os.path.abspath(input_pdb)
        input_fasta_abs = os.path.abspath(input_fasta)
    
        os.chdir(working_dir)
    
        # Use PSIPRED to generate secondary structure predictions
   
        ss2_file = self.run_psipred(input_fasta_abs)
   
        # Use CS-BLAST to generate the sequence profile

        checkpoint_binary = self.run_csblast(input_fasta_abs)
        checkpoint_text = checkpoint_binary[:-5] + 'checkpoint' 
        
        sequence = IO.fasta_file_to_sequences(input_fasta_abs)[1][0]
        self.parse_checkpoint_file(checkpoint_binary, checkpoint_text, len(sequence))

        # Create the PSSM file
        
        pssm_file = self.create_pssm(input_fasta_abs)
        
        # Get the phipsi file

        phipsi_file = self.run_sparksx(input_fasta_abs, pssm_file)
        
        # Do fragment picking
   
        pick_cmd = [self.fragment_picker_cmd]
        
        if self.rosetta_database:
            pick_cmd += ['-in:path:database', self.rosetta_database]
        
        # Pick 3-mer fragment library:
        pick_cmd += ['-in::file::vall', self.vall_abs,
                    '-in::file::fasta', input_fasta_abs,
                    '-in::file::s', input_pdb_abs,
                    '-frags::scoring::config', self.weights_file_abs,
                    '-frags::frag_sizes', '3', 
                    '-frags::n_candidates', '1000',
                    '-frags::n_frags', '200',
                    '-frags::ss_pred', ss2_file, 'predA',
                    '-in:file:checkpoint', checkpoint_text,
                    '-spine_x', phipsi_file,
                    '-out::file::frag_prefix', 'frags',
                    '-frags::describe_fragments', 'frags.fsc']
        
        subprocess.check_call(pick_cmd)
        
        # Pick 9-mer fragment library:
        pick_cmd.remove('-frags::frag_sizes')
        pick_cmd.remove('3')
        pick_cmd += ['-frags::frag_sizes', '9']
        
        subprocess.check_call(pick_cmd)
        
        # Clean the temporary files

        self.clean(input_fasta)

        # Go back to the previous working directory
    
        os.chdir(cwd)