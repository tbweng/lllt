import os, glob
from nipype import Node, Workflow  # components to construct workflow
import nipype.interfaces.fsl as fsl # importing FSL interface functions
from nipype.interfaces.io import DataSink  # datasink

dataDir = '/scratch/05201/bhickson/preprocdata'
indCope = '4'  # the contrast of interest, 4: 2B > others, 6: 2B > 0B
nreg_id = 'noGSR'
task_id = 'nback'

###########
#
# SETTING UP THE SECOND LEVEL ANALYSIS NODES (combine runs within a subject with a fixed-effects model)
#
###########

###########
#
# A LIST OF COPE, VARCOPE, AND MASK FILES TO BE MEREGED
#
###########
# directory where level 1 is
baseDir = os.path.join(dataDir, 'feat_dir')


for ses in ['A', 'B']:
    outDir = os.path.join(dataDir, 'level2', nreg_id, ses, 'cope', indCope)

    # a list of subjects
    subject_list = [ i.split('_')[3].split('-')[1] for i in sorted(glob.glob(baseDir + '/' + nreg_id + '/run-01/ses-_' + ses + '_sub*')) ]

    listCopeFiles = []
    listVarcopeFiles = []
    listMaskFiles = []

    for iSubj in subject_list:
        for run_id in ['01', '02']:
            
            # full path to a cope image
            pathCope = os.path.join(baseDir,
                                    nreg_id,
                                    ('run-' + run_id +
                                    '/ses-_' + ses +
                                    '_sub-' + iSubj +
                                    '_task-' + task_id),
                                    'run0.feat',
                                    'stats',
                                    'cope' + indCope + '.nii.gz')
            listCopeFiles.append(pathCope)

            # full path to a varcope image
            pathVarcope = os.path.join(baseDir,
                                    nreg_id,
                                    ('run-' + run_id +
                                        '/ses-_' + ses +
                                        '_sub-' + iSubj +
                                        '_task-' + task_id),
                                    'run0.feat',
                                    'stats',
                                    'varcope' + indCope + '.nii.gz')
            listVarcopeFiles.append(pathVarcope)

            # full path to a mask image
            pathMask = os.path.join(baseDir,
                                    nreg_id,
                                    ('run-' + run_id +
                                    '/ses-_' + ses +
                                    '_sub-' + iSubj +
                                    '_task-' + task_id),
                                    'run0.feat',
                                    'mask.nii.gz')
            listMaskFiles.append(pathMask)

    ###########
    #
    # SETTING UP THE SECOND LEVEL ANALYSIS NODES
    #
    ###########

    # Dictionary with regressors
    vector_length = len(subject_list)*2

    evs = vector_length//2

    dictReg = {
        'ev%03d' % i : ([0] * (i * 2) + [1, 1] + [0] * vector_length)[:vector_length]
        for i in range(evs)
    }

    # Contrasts (list of lists)
    contrastList = []

    contrastDict = {
        'cont%03d' % i: ([0] * i + [1] + [0] * evs)[:evs]
        for i in range(evs)
    }

    for key, value in contrastDict.items():
        temp = [key, 'T', list(dictReg.keys()), value]
        contrastList.append(temp)

    # Setting up the second level analysis model node
    level2design = Node(fsl.MultipleRegressDesign(contrasts=contrastList,
                                                regressors=dictReg),
                        name='level2design')

    # Model calculation by FLAMEO
    flameo = Node(fsl.FLAMEO(run_mode='fe'),
                name="flameo")



    ###########
    #
    # NODES FOR THE MERGING IMAGES
    #
    ###########
    # merging cope files
    copemerge = Node(fsl.Merge(dimension='t',
                            in_files=listCopeFiles),
                    name="copemerge")

    # merging varcope files
    varcopemerge = Node(fsl.Merge(dimension='t',
                            in_files=listVarcopeFiles),
                        name="varcopemerge")

    # merging mask files
    maskmerge = Node(fsl.Merge(dimension='t',
                            in_files=listMaskFiles),
                    name="maskmerge")

    # calculating the minimum across time points on merged mask image
    minmask = Node(fsl.MinImage(),
                name="minmask")


    # creating datasink to collect outputs
    datasink = Node(DataSink(base_directory=
                            os.path.join(outDir, 'L2')),
                    name='datasink')



    ###########
    #
    # SETTING UP THE WORKFLOW NODES
    #
    ###########

    # creating the workflow
    secondLevel = Workflow(name="Level2_wf", base_dir=outDir)

    # connecting nodes
    secondLevel.connect(level2design, 'design_mat', flameo, 'design_file')
    secondLevel.connect(level2design, 'design_con', flameo, 't_con_file')
    secondLevel.connect(level2design, 'design_grp', flameo, 'cov_split_file')
    secondLevel.connect(copemerge, 'merged_file', flameo, 'cope_file')
    secondLevel.connect(varcopemerge, 'merged_file', flameo, 'var_cope_file')
    secondLevel.connect(maskmerge, 'merged_file', minmask, 'in_file')
    secondLevel.connect(minmask, 'out_file', flameo, 'mask_file')
    secondLevel.connect(flameo, 'stats_dir', datasink, 'stats_dir')

    # write out graph
    # secondLevel.write_graph(graph2use='orig', dotfilename='graph_orig.dot')

    # running the workflow
    secondLevel.run()
