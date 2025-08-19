import React, { createContext, useContext, useState, useCallback, useEffect } from "react";

// Define the context type
export type EnvironmentContextType = {
  selectedEnv: string | null;
  setSelectedEnv: (env: string) => void;
  environments: any[];
  setEnvironments: React.Dispatch<React.SetStateAction<any[]>>;
  loadingEnvironments: boolean;
  setLoadingEnvironments: React.Dispatch<React.SetStateAction<boolean>>;
  refetchEnvironments: () => Promise<any[]>;
  refreshEnvFiles: (envId: string) => Promise<any[]>;
};

// Create the context
const EnvironmentContext = createContext<EnvironmentContextType | undefined>(undefined);

// Provider component
export const EnvironmentProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log('[EnvironmentProvider] EnvironmentProvider mounting');
  const [selectedEnv, setSelectedEnv] = useState<string | null>(null);
  const [environments, setEnvironments] = useState<any[]>([]);
  const [loadingEnvironments, setLoadingEnvironments] = useState<boolean>(false);

  useEffect(() => {
    console.log('[EnvironmentProvider] EnvironmentProvider mounted');
    return () => {
      console.log('[EnvironmentProvider] EnvironmentProvider unmounting');
    };
  }, []);

  const refetchEnvironments = useCallback(() => {
    console.log('[EnvironmentProvider] Refetching environments...');
    setLoadingEnvironments(true);
    return fetch('/api/environments', { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        console.log('[EnvironmentProvider] Environments fetched:', data.environments);
        setEnvironments(data.environments || []);
        return data.environments || [];
      })
      .finally(() => setLoadingEnvironments(false));
  }, []);

  const refreshEnvFiles = useCallback((envId: string) => {
    console.log('[EnvironmentProvider] Refreshing files for environment:', envId);
    return fetch(`/api/environments/${envId}/files`, {
      credentials: 'include',
      cache: 'no-store'
    })
      .then(res => res.json())
      .then(data => {
        console.log('[EnvironmentProvider] Files refreshed for environment:', envId, data.files);
        return data.files || [];
      });
  }, []);

  console.log('[EnvironmentProvider] EnvironmentProvider rendering with state:', {
    selectedEnv,
    environmentsCount: environments.length,
    loadingEnvironments
  });

  return (
    <EnvironmentContext.Provider value={{ selectedEnv, setSelectedEnv, environments, setEnvironments, loadingEnvironments, setLoadingEnvironments, refetchEnvironments, refreshEnvFiles }}>
      {children}
    </EnvironmentContext.Provider>
  );
};

// Hook to use the context
export const useEnvironment = () => {
  const context = useContext(EnvironmentContext);
  if (!context) {
    throw new Error("useEnvironment must be used within an EnvironmentProvider");
  }
  return context;
}; 