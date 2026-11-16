/* BlueMolecule — drop-down menus (issue #2).
 *
 * The system pop-up of a <select> does not stay open after a click in some windows (for example the Qt window on
 * KDE Wayland): one has to press, drag to the option and release. So BlueMolecule opens its own menu for every
 * single-choice <select> instead. The <select> elements stay as they are: they show the current value, keep their
 * styles and change events, and code (or a test) can still set their value directly.
 *
 *   - mouse: click to open, click an option; or press, drag to an option and release, as with a system menu;
 *   - keyboard on a focused drop-down: Enter, Space, F4 or Alt+Down open it; in the menu ↑ ↓ Home End PageUp
 *     PageDown move, Enter or Space choose, Esc or Tab close, letters jump to the next option starting with them;
 *   - the menu closes on a click elsewhere, on scrolling outside it, on resize and when the window loses focus;
 *   - after a mouse choice the drop-down gives the keyboard back to the window (outside dialogs), so the arrow keys
 *     step through the geometries again.
 * A <select data-native> keeps the system pop-up.
 */
(function () {
  'use strict';
  const G = typeof window !== 'undefined' ? window : globalThis;
  const BM = (G.BM = G.BM || {});
  if (typeof document === 'undefined') { BM.selectMenu = { open() {}, close() {}, isOpen: () => false }; return; }

  let cur = null;   // { sel, menu, items, hi, armed, x, y, t, typed, typedAt }

  function usable(el) {
    return el && el.tagName === 'SELECT' && !el.multiple && !(el.size > 1) && !el.disabled && !el.hasAttribute('data-native');
  }

  function build(sel) {
    const menu = document.createElement('div');
    menu.className = 'bm-menu';
    menu.setAttribute('role', 'listbox');
    const name = sel.getAttribute('aria-label') || sel.title;
    if (name) menu.setAttribute('aria-label', name);
    const items = [];
    const add = (o, grouped, groupOff) => {
      const el = document.createElement('div');
      el.className = 'bm-opt' + (grouped ? ' grouped' : '') + (o.disabled || groupOff ? ' off' : '');
      el.setAttribute('role', 'option');
      el.textContent = o.label || o.textContent;
      if (o.title) el.title = o.title;
      if (o.selected) el.setAttribute('aria-selected', 'true');
      el._opt = o;
      el._off = !!(o.disabled || groupOff);
      items.push(el);
      menu.appendChild(el);
    };
    for (const ch of sel.children) {
      if (ch.tagName === 'OPTGROUP') {
        const h = document.createElement('div');
        h.className = 'bm-grp';
        h.textContent = ch.label;
        menu.appendChild(h);
        for (const o of ch.children) if (o.tagName === 'OPTION') add(o, true, ch.disabled);
      } else if (ch.tagName === 'OPTION') add(ch, false, false);
    }
    return { menu, items };
  }

  function place(sel, menu) {
    const r = sel.getBoundingClientRect(), vw = window.innerWidth, vh = window.innerHeight;
    menu.style.minWidth = Math.max(120, r.width) + 'px';
    menu.style.maxHeight = '';
    const below = vh - r.bottom - 8, above = r.top - 8;
    const h = Math.min(menu.scrollHeight, 360);
    const up = h > below && above > below;
    const room = Math.max(80, up ? above : below);
    menu.style.maxHeight = Math.min(360, room) + 'px';
    const mh = Math.min(h, room);
    menu.style.top = (up ? r.top - mh - 2 : r.bottom + 2) + 'px';
    const w = menu.offsetWidth;
    menu.style.left = Math.max(4, Math.min(r.left, vw - w - 4)) + 'px';
  }

  function highlight(k, scroll) {
    if (!cur) return;
    if (cur.hi >= 0 && cur.items[cur.hi]) cur.items[cur.hi].classList.remove('hi');
    cur.hi = k;
    const el = cur.items[k];
    if (!el) return;
    el.classList.add('hi');
    cur.menu.setAttribute('aria-activedescendant', el.id || (el.id = 'bm-opt-' + k));
    if (scroll) {   // scroll the menu only (never the panels behind it)
      const m = cur.menu;
      if (el.offsetTop < m.scrollTop) m.scrollTop = el.offsetTop - 4;
      else if (el.offsetTop + el.offsetHeight > m.scrollTop + m.clientHeight) m.scrollTop = el.offsetTop + el.offsetHeight - m.clientHeight + 4;
    }
  }

  function step(from, dir) {
    const n = cur.items.length;
    for (let k = from + dir; k >= 0 && k < n; k += dir) if (!cur.items[k]._off) return k;
    return from;
  }

  function open(sel, how) {
    if (!usable(sel)) return;
    close();
    const { menu, items } = build(sel);
    if (!items.length) return;
    // inside a modal dialog (drawn above the whole page) the menu must be part of the dialog to be seen; being
    // position: fixed it is still placed against the window and not cut by the dialog's edges
    (sel.closest('dialog[open]') || document.body).appendChild(menu);
    cur = { sel, menu, items, hi: -1, armed: how !== 'mouse', x: 0, y: 0, t: Date.now(), typed: '', typedAt: 0 };
    sel.classList.add('bm-open');
    sel.setAttribute('aria-expanded', 'true');
    place(sel, menu);
    const k = items.findIndex(el => el._opt.selected);
    highlight(k >= 0 ? k : step(-1, 1), true);
    menu.addEventListener('mousedown', e => e.preventDefault());   // keep the focus on the <select>
  }

  function close() {
    if (!cur) return;
    const c = cur;
    cur = null;
    c.menu.remove();
    c.sel.classList.remove('bm-open');
    c.sel.setAttribute('aria-expanded', 'false');
  }

  function choose(k, how) {
    if (!cur) return;
    const el = cur.items[k], sel = cur.sel;
    if (!el || el._off) return;
    close();
    // after a choice made with the mouse, the keyboard goes back to the window (arrow keys step through the
    // geometries, letters are shortcuts) unless the drop-down is in a dialog; after a keyboard choice it stays here
    if (how === 'mouse' && !sel.closest('dialog')) sel.blur();
    else { try { sel.focus({ preventScroll: true }); } catch (e) { sel.focus(); } }
    if (sel.value === el._opt.value && el._opt.selected) return;
    el._opt.selected = true;
    sel.dispatchEvent(new Event('input', { bubbles: true }));
    sel.dispatchEvent(new Event('change', { bubbles: true }));
  }

  // ---- mouse: open instead of the system pop-up; click or drag-and-release to choose
  document.addEventListener('mousedown', e => {
    if (cur && cur.menu.contains(e.target)) return;
    const sel = e.target && e.target.closest ? e.target.closest('select') : null;
    if (cur && cur.sel === sel) { e.preventDefault(); close(); return; }
    close();
    if (!usable(sel) || e.button !== 0) return;
    e.preventDefault();   // no system pop-up
    try { sel.focus({ preventScroll: true }); } catch (err) { sel.focus(); }
    open(sel, 'mouse');
    if (cur) { cur.x = e.clientX; cur.y = e.clientY; cur.t = Date.now(); }
  }, true);
  // while the button is held after pressing on the <select>, the browser sends the mouse events to the <select>,
  // so the option under the pointer is found by position
  const optionAt = e => {
    const under = document.elementFromPoint(e.clientX, e.clientY);
    const el = under && under.closest ? under.closest('.bm-opt') : null;
    return el && cur && cur.menu.contains(el) ? el : null;
  };
  document.addEventListener('mousemove', e => {
    if (!cur) return;
    const el = optionAt(e);
    if (el && !el._off) highlight(cur.items.indexOf(el), false);
  }, true);
  document.addEventListener('mouseup', e => {
    if (!cur || e.button !== 0) return;
    const el = optionAt(e);
    const dragged = Math.abs(e.clientX - cur.x) + Math.abs(e.clientY - cur.y) > 4 || Date.now() - cur.t > 350;
    if (el && cur.menu.contains(el) && (cur.armed || dragged)) { choose(cur.items.indexOf(el), 'mouse'); return; }
    cur.armed = true;   // the click that opened the menu: the next click chooses
  }, true);

  // ---- keyboard
  document.addEventListener('keydown', e => {
    if (cur) {
      const n = cur.items.length;
      let k = null;
      switch (e.key) {
        case 'ArrowDown': k = step(cur.hi, 1); break;
        case 'ArrowUp': k = e.altKey ? null : step(cur.hi, -1); if (e.altKey) { e.preventDefault(); e.stopPropagation(); close(); return; } break;
        case 'Home': k = step(-1, 1); break;
        case 'End': k = step(n, -1); break;
        case 'PageDown': k = Math.min(n - 1, cur.hi + 10); if (cur.items[k]._off) k = step(k, -1); break;
        case 'PageUp': k = Math.max(0, cur.hi - 10); if (cur.items[k]._off) k = step(k, 1); break;
        case 'Enter': case ' ':
          e.preventDefault(); e.stopPropagation(); choose(cur.hi, 'key'); return;
        case 'Escape':
          e.preventDefault(); e.stopPropagation(); close(); return;
        case 'Tab': close(); return;
        default:
          if (e.key.length === 1 && !e.ctrlKey && !e.metaKey && !e.altKey) {   // type to jump
            const now = Date.now();
            cur.typed = (now - cur.typedAt < 700 ? cur.typed : '') + e.key.toLowerCase();
            cur.typedAt = now;
            const at = (j) => cur.items[j].textContent.trim().toLowerCase();
            for (let d = cur.typed.length > 1 ? 0 : 1; d <= n; d++) {
              const j = (cur.hi + d) % n;
              if (!cur.items[j]._off && at(j).startsWith(cur.typed)) { k = j; break; }
            }
          }
      }
      if (k !== null) { e.preventDefault(); e.stopPropagation(); highlight(k, true); }
      else if (e.key.length === 1) { e.preventDefault(); e.stopPropagation(); }
      return;
    }
    const sel = e.target;
    if (!usable(sel)) return;
    if (e.key === ' ' || e.key === 'Enter' || e.key === 'F4' || (e.altKey && (e.key === 'ArrowDown' || e.key === 'ArrowUp'))) {
      e.preventDefault();
      e.stopPropagation();
      open(sel, 'key');
    }
  }, true);

  // ---- close when the page moves under the menu or the window loses focus
  document.addEventListener('scroll', e => { if (cur && !cur.menu.contains(e.target)) close(); }, true);
  window.addEventListener('resize', () => close());
  window.addEventListener('blur', () => close());
  document.addEventListener('focusin', e => { if (cur && e.target !== cur.sel && !cur.menu.contains(e.target)) close(); });

  BM.selectMenu = { open: sel => open(sel, 'key'), close, isOpen: () => !!cur, current: () => (cur ? cur.sel : null) };
})();
