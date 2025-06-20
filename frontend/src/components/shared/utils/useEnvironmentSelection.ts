import { useEffect } from 'react';
import { useSearchParams } from 'react-router-dom';

interface Environment {
  environment_id: string;
  title: string;
}

export const useEnvironmentSelection = (
  environments: Environment[],
  selectedEnv?: string,
  setSelectedEnv?: (envId: string) => void
) => {
  const [searchParams, setSearchParams] = useSearchParams();

  // Load environment from URL
  useEffect(() => {
    if (!setSelectedEnv) return;
    const envId = searchParams.get('env');
    if (envId && environments.some(e => e.environment_id === envId)) {
      setSelectedEnv(envId);
    } else if (environments.length > 0) {
      setSelectedEnv(environments[0].environment_id);
      setSearchParams({ env: environments[0].environment_id }, { replace: true });
    }
  }, [environments, setSelectedEnv, searchParams]);

  // Update URL when environment changes
  useEffect(() => {
    if (selectedEnv) {
      setSearchParams({ env: selectedEnv }, { replace: true });
    }
  }, [selectedEnv, setSearchParams]);

  return { selectedEnv };
}; 