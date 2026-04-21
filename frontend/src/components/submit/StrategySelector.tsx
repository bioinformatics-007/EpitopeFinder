"use client";

import React, { useEffect, useState } from 'react';
import { useWizardStore } from '@/lib/store';
import { Tooltip } from '@/components/ui/Tooltip';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/Card';
import { Badge } from '@/components/ui/Badge';
import { api } from '@/lib/api';
import type { StrategyInfo } from '@/lib/types';

export function StrategySelector() {
  const { strategy, setStrategy } = useWizardStore();
  const [strategies, setStrategies] = useState<StrategyInfo[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    api.getStrategies()
      .then(setStrategies)
      .catch(console.error)
      .finally(() => setLoading(false));
  }, []);

  return (
    <div className="space-y-8">
      <div className="flex flex-col space-y-3">
        <h2 className="text-2xl font-bold text-stone-900 ">Select Strategy</h2>
        <p className="text-stone-500 font-medium leading-relaxed">
          Choose an analysis pipeline strategy. Our platform offers 6 specialized biological assessment 
          pathways tailored for different research objectives.
        </p>
      </div>

      {loading ? (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {[1, 2, 3, 4, 5, 6].map((i) => (
            <div key={i} className="h-48 rounded-sm skeleton" />
          ))}
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {strategies.map((info) => {
            const isSelected = strategy === info.number;
            
            return (
              <Card
                key={info.number}
                hoverable
                onClick={() => setStrategy(info.number)}
                className={`relative overflow-hidden group ${
                  isSelected 
                    ? 'border-amber-600 bg-amber-50/30 ' 
                    : 'bg-white'
                }`}
              >
                {isSelected && (
                  <div className="absolute top-0 right-0 p-5">
                    <div className="h-3 w-3 rounded-sm bg-amber-600 0_0_12px_rgba(16,185,129,0.8)] " />
                  </div>
                )}
                
                <CardHeader className="pb-3 border-b border-transparent group-">
                  <div className="flex items-center gap-4">
                    <div className={`flex h-10 w-10 items-center justify-center rounded-sm font-mono text-base font-bold border-2 ${
                      isSelected 
                        ? 'bg-amber-600 text-white border-amber-600 rotate-12 scale-110' 
                        : 'bg-stone-50 text-stone-400 border-stone-100 group-'
                    }`}>
                      {info.number}
                    </div>
                    <CardTitle className={`flex items-center text-lg ${isSelected ? 'text-amber-900' : 'text-stone-900'}`}>
                      {info.name}
                      {info.number === 1 && <Tooltip text="Comprehensive prediction of B-cell, CTL, and HTL epitopes from a target sequence." />}
                      {info.number === 2 && <Tooltip text="Focuses on a specific type of epitope (B-cell, CTL, HTL) or a specific pathogen." />}
                      {info.number === 3 && <Tooltip text="Specialized pipeline for Toxins, Allergens, or Gut Flora compatibility." />}
                      {info.number === 4 && <Tooltip text="Evaluate already predicted epitopes across multiple criteria (Toxicity, Allergenicity, etc.)." />}
                      {info.number === 5 && <Tooltip text="Build vaccine sequence using your previously generated epitope outputs instead of predicting from scratch." />}
                      {info.number === 6 && <Tooltip text="Perform structural modeling, assemble epitopes, and analyze the final vaccine construct." />}
                    </CardTitle>
                  </div>
                </CardHeader>
                <CardContent className="pt-4">
                  <CardDescription className={`leading-relaxed mb-6 font-medium ${isSelected ? 'text-amber-800/80' : 'text-stone-500'}`}>
                    {info.description}
                  </CardDescription>
                  <div className="flex flex-wrap gap-2">
                    {info.available_tools.slice(0, 4).map((tool) => (
                      <Badge key={tool} variant={isSelected ? 'success' : 'default'} className="px-2 py-0.5">
                        {tool}
                      </Badge>
                    ))}
                    {info.available_tools.length > 4 && (
                      <Badge variant={isSelected ? 'success' : 'default'}>
                        +{info.available_tools.length - 4} more
                      </Badge>
                    )}
                  </div>
                </CardContent>
              </Card>
            );
          })}
        </div>
      )}
    </div>
  );
}
