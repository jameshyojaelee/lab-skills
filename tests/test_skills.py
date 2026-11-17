"""Tests on synthetic data with known answers. Run from the repo root:
    python3 tests/test_skills.py        (or: python3 -m pytest -q tests)"""
import os, sys, tempfile
import numpy as np

os.environ["LABSTYLE_STRICT"] = "1"   # any figure-rule violation fails the test

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in ('skills/gel-quant/scripts', 'skills/qpcr-quant/scripts', 'skills/tapestation-qc/scripts',
          'skills/cell-count/scripts', 'skills/library-quant/scripts', 'tests'):
    sys.path.insert(0, os.path.join(ROOT, p))
import gel_quant, qpcr_quant, tapestation_qc, cell_count, library_quant
import synthetic_gel, synthetic_qpcr, synthetic_tapestation, synthetic_countess, synthetic_cfx


def test_gel_synthetic():
    cfg = os.path.join(ROOT, 'skills', 'gel-quant', 'examples', 'synthetic', 'config.yaml')
    with tempfile.TemporaryDirectory() as td:
        img = os.path.join(td, 'gel.tif')
        synthetic_gel.render(img)
        g, summ = gel_quant.run(cfg, image=img, outdir=td, plots=True)
        assert os.path.exists(os.path.join(td, 'metrics.png'))
    for r in g.metric_rows:
        want = synthetic_gel.truth(r['lane'])[r['metric']]
        assert abs(r['value'] - want) < 0.02, (r, want)
    st = summ['standards']
    assert st['usable'] == 1 and any(p['clip'] for p in st['points'])   # high standard saturates


def test_gel_committed_example_matches_generator():
    """The committed example image must be what the generator makes (deterministic)."""
    from PIL import Image
    path = os.path.join(ROOT, 'skills', 'gel-quant', 'examples', 'synthetic', 'synthetic_gel.tif')
    assert np.array_equal(np.array(Image.open(path)), synthetic_gel.render())


def test_qpcr_synthetic():
    cfg = os.path.join(ROOT, 'skills', 'qpcr-quant', 'examples', 'single-bead-bulk', 'config.yaml')
    with tempfile.TemporaryDirectory() as td:
        xl = os.path.join(td, 'synthetic.xlsx')
        truth = synthetic_qpcr.make(xl)
        s = qpcr_quant.run(cfg, export=xl, outdir=td, make_plots=True)
        assert os.path.exists(os.path.join(td, 'qpcr_qc.png'))
    f = s['standard_curve']
    assert abs(f['slope'] - truth['slope']) < 0.05 and f['r2'] > 0.999
    sb = s['single_bead']['0.2_per_well']
    assert abs(np.log10(sb['per_bead_median'] / truth['unit'])) < 0.1
    assert abs(np.log10(sb['background_copies'] / truth['bg'])) < 0.1
    assert sb['bead_wells'] == int((truth['k'] > 0).sum())
    for k, v in truth['bulk'].items():
        assert abs(np.log10(s['bulk'][k]['mean'] / v)) < 0.02, (k, s['bulk'][k], v)


def test_auto_threshold():
    ct = np.array([20.1, 20.0, 19.9, 20.2, 16.4, 16.6, 16.5])
    assert 16.6 < qpcr_quant.auto_threshold(ct) < 19.9


def test_tapestation_synthetic():
    with tempfile.TemporaryDirectory() as td:
        truth = synthetic_tapestation.write(os.path.join(td, 'run'))
        o = tapestation_qc.run(os.path.join(td, 'run'), outdir=os.path.join(td, 'out'))
        assert os.path.exists(os.path.join(td, 'out', 'gel_labeled.png'))
    by = {r['well']: r for r in o['samples']}
    for w, (bp, pa, dbp, dpa, rt, mol) in truth.items():
        assert by[w]['main_peak_bp'] == bp and abs(by[w]['region_molarity_nM'] - mol / 1000) < 1e-9
        assert abs(by[w]['pct_area_below'] - (dpa if dbp else 0)) < 1e-9
        assert abs(o['size_check'][w]['trace_bp'] - bp) / bp < 0.03      # ladder sizing from the trace
    assert any('B1' in x and 'below 200' in x for x in o['warnings'])
    assert any('C1' in x and 'only 41' in x for x in o['warnings'])


def test_countess_synthetic():
    with tempfile.TemporaryDirectory() as td:
        truth = synthetic_countess.write(os.path.join(td, 'cc'))
        recs, groups, warn = cell_count.run([os.path.join(td, 'cc')], seed_cells=2e5, outdir=os.path.join(td, 'out'))
    assert len(recs) == len(truth)
    for (name, n, live, k), r in zip(truth, recs):
        assert r['sample'] == name and r['json_pos'] == k and r['cube1_counted'] == k
        assert abs(r['viability_pct'] - 100 * round(n * live) / n) < 0.01
        assert r['cube1_ci_lo'] < 100 * k / n < r['cube1_ci_hi']
    g = {x['group']: x for x in groups}
    assert g['cellsA']['n_reads'] == 2 and abs(g['cellsA']['cube1_pooled_pct'] - 100 * 230 / 780) < 1e-6
    assert abs(g['cellsB']['ul_for_seed'] - 2e5 / (400 * 2500) * 1000) < 1e-6


def test_library_quant_synthetic():
    cfg = os.path.join(ROOT, 'skills', 'library-quant', 'examples', 'synthetic', 'config.yaml')
    with tempfile.TemporaryDirectory() as td:
        truth = synthetic_cfx.write(os.path.join(td, 'export'))
        fit, libs, warn = library_quant.run(cfg, export_dir=os.path.join(td, 'export'), outdir=td)
        assert os.path.exists(os.path.join(td, 'library_quant.png'))
    assert abs(fit['slope'] - synthetic_cfx.SLOPE) < 0.01 and fit['r2'] > 0.999
    for r in libs:
        nM, bp = truth[r['name']]
        assert abs(r['nM_size_corrected'] / nM - 1) < 0.05, (r['name'], r['nM_size_corrected'], nM)
        assert abs(r['delta_cq'] - r['delta_cq_expected']) < 0.2
    l1 = [r for r in libs if r['name'] == 'L1'][0]
    assert abs(l1['qpcr_over_fluorometric'] - 1) < 0.05


def test_countess_timecourse_config():
    import json
    import cell_count_timecourse
    with tempfile.TemporaryDirectory() as td:
        truth = synthetic_countess.write(os.path.join(td, 'raw'))
        import glob
        for source in glob.glob(os.path.join(td, 'raw', '*', '*.csv')):
            with open(source) as fh:
                content = fh.read()
            with open(source, 'w') as fh:
                fh.write(content.replace(':00,', ':00.1234567,', 1))
        cfg = dict(input_globs=['raw/*'], arm_regex=r'(?P<arm>cells[A-B])',
                   arm_labels={'cellsA': 'Arm A', 'cellsB': 'Arm B'}, channel='GFP',
                   start_time='2026-01-01T00:00:00', title='Synthetic time course',
                   repeat_count_ids=['1001'], outdir='results')
        config = os.path.join(td, 'config.json')
        with open(config, 'w') as fh:
            json.dump(cfg, fh)
        rows, warnings = cell_count_timecourse.run(config)
        assert len(rows) == len(truth)
        for row, (_, n, _, k) in zip(rows, truth):
            assert row['objects'] == n and row['positive'] == k
            assert abs(row['positive_pct'] - 100*k/n) < 1e-10
            assert row['nonclumped_objects'] == n
            assert row['ci_low'] < row['positive_pct'] < row['ci_high']
        assert rows[1]['read_role'] == 'repeat'
        assert abs(rows[0]['hours_post_start'] - (10 + 0.123456 / 3600)) < 1e-10
        assert os.path.exists(os.path.join(td, 'results', 'timecourse.png'))
        # Corrupt a second export's ID: duplicate exports must never inflate the series.
        source = rows[1]['source_csv']
        with open(source) as fh:
            content = fh.read()
        with open(source, 'w') as fh:
            fh.write(content.replace('1001,', '1000,', 1))
        try:
            cell_count_timecourse.run(config, plots=False)
        except ValueError as exc:
            assert 'Duplicate Count ID' in str(exc)
        else:
            raise AssertionError('Duplicate ID was accepted')


def test_default_font_is_helvetica():
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(ROOT, '_shared'))
    import labstyle
    labstyle.setup()
    order = plt.rcParams['font.sans-serif']
    assert order[0] == 'Helvetica'
    # Arial and Liberation Sans must come after every Helvetica-family name.
    last_helv = max(i for i, n in enumerate(order) if n in labstyle.HELVETICA_FAMILY)
    assert last_helv < order.index('Arial') < order.index('Liberation Sans')
    # The PDF must embed the font that was resolved (editable text, Type 42).
    name = labstyle.resolved_font()
    with tempfile.TemporaryDirectory() as td:
        fig, ax = labstyle.figure('single')
        ax.plot([0, 1], [0, 1]); ax.set_xlabel('x (unit)'); ax.set_ylabel('y (unit)')
        labstyle.save(fig, os.path.join(td, 'f'))
        raw = open(os.path.join(td, 'f.pdf'), 'rb').read()
        assert name.replace(' ', '').encode() in raw, (name, 'not embedded in the PDF')


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print('ok', name)
