import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework

process = cepgen.Module('grape',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        LPAIR = PDG.electron,
    ),
    inKinematics = cepgen.Parameters(
        pz = (27.52, 820.),
        pdgIds = (11, 2212),
    ),
)

generator = cepgen.Parameters(
    numEvents = 100
)
text = cepgen.Module('text',
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 250, 10)]),
        'm(ob1)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
output = cepgen.Sequence(text)

