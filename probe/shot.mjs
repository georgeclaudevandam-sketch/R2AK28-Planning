import { chromium } from 'playwright';
const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
await p.goto('http://localhost:8765/index.html'); await p.waitForTimeout(6000);
for (const [n, la, lo, z] of [['nanaimo11', 49.2, -123.8, 11], ['nanaimo12', 49.15, -123.85, 12], ['gulf11', 48.85, -123.35, 11], ['vic12', 48.42, -123.37, 12]]) {
  await p.evaluate(([la, lo, z]) => { setFollow(false); map.setView([la, lo], z); }, [la, lo, z]);
  await p.waitForTimeout(12000);
  await p.screenshot({ path: `probe/${n}.png`, clip: { x: 0, y: 0, width: 760, height: 900 } });
}
await b.close();
