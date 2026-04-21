"use client";

import React from 'react';
import Link from 'next/link';
import { usePathname } from 'next/navigation';

export function HorizontalMenu() {
  const pathname = usePathname();

  const navItems = [
    { name: 'Home', href: '/' },
    { name: 'Design Vaccine', href: '/submit' },
    { name: 'Results Viewer', href: '/results' },
    { name: 'Algorithm', href: '/algorithm' },
    { name: 'Help', href: '/help' },
    { name: 'Contact', href: '/contact' },
  ];

  return (
    <nav className="w-full bg-[#635B53] border-b border-[#3A332D]">
      <div className="flex flex-wrap">
        {navItems.map((item) => {
          const isActive = item.href === '/' 
            ? pathname === '/' 
            : pathname.startsWith(item.href);
            
          return (
            <Link
              key={item.name}
              href={item.href}
              className={`px-4 py-1.5 text-[13px] font-bold border-r border-[#4A433A] ${
                isActive 
                  ? 'bg-[#3A332D] text-[#ECEF01]' 
                  : 'text-white hover:text-[#ECEF01] hover:bg-[#726A61]'
              } `}
            >
              {item.name}
            </Link>
          );
        })}
      </div>
    </nav>
  );
}
