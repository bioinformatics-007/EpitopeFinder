import React from 'react';
import Link from 'next/link';

export function Footer() {
  return (
    <footer className="border-t border-[#999] bg-[#3A332D] text-white text-center py-3 text-sm">
      <p>
        <Link href="/algorithm" className="hover:text-[#ECEF01]">Algorithm</Link>
        {' | '}
        <Link href="/help" className="hover:text-[#ECEF01]">Help</Link>
        {' | '}
        <Link href="/contact" className="hover:text-[#ECEF01]">Contact Us</Link>
        {' | '}
        <a href="https://github.com/yuktika12/vaccine-design-pipeline" target="_blank" rel="noreferrer" className="hover:text-[#ECEF01]">GitHub</a>
      </p>
    </footer>
  );
}
