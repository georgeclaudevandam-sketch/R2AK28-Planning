import { chromium } from 'playwright';
const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
const views = [['gulf', 48.85, -123.35, 11], ['gulf10', 48.85, -123.35, 10], ['vic12', 48.42, -123.37, 12]];
for (const [dpi, ts] of [[192, 1.3], [96, 2], [64, 2.6]]) {
  await p.goto(`http://localhost:8765/index.html?chartdpi=${dpi}&charttext=${ts}`); await p.waitForTimeout(6000);
  for (const [n, la, lo, z] of views) {
    await p.evaluate(([la, lo, z]) => { setFollow(false); map.setView([la, lo], z); }, [la, lo, z]);
    await p.waitForTimeout(10000);
    await p.screenshot({ path: `probe/${n}-dpi${dpi}.png`, clip: { x: 0, y: 0, width: 760, height: 900 } });
  }
}
await b.close();
