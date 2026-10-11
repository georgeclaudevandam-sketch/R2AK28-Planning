import { chromium } from 'playwright';
const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
p.on('console', m => { if (m.type()==='error') console.log('console', m.text()); });
await p.goto('http://localhost:8765/index.html'); await p.waitForTimeout(6000);
const views = [['victoria', 48.42, -123.37, 12], ['gulf', 48.85, -123.35, 11], ['haro', 48.6, -123.2, 11], ['porlier', 49.01, -123.58, 12]];
for (const [n, la, lo, z] of views) {
  for (const mode of ['default', 'chs-only']) {
    await p.evaluate(([la, lo, z, mode]) => { setFollow(false); map.eachLayer(l => { if (l._noaaTag && mode === 'chs-only') map.removeLayer(l); }); map.setView([la, lo], z); }, [la, lo, z, mode]);
    await p.waitForTimeout(9000);
    await p.screenshot({ path: `probe/${n}-${mode}.png`, clip: { x: 0, y: 0, width: 1000, height: 900 } });
  }
  await p.reload(); await p.waitForTimeout(6000);
}
await b.close();
// rerun
