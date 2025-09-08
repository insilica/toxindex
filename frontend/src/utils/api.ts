/**
 * Utility functions for making authenticated API requests
 */

/**
 * Make an authenticated API request
 */
export const authenticatedFetch = async (
  url: string,
  options: RequestInit = {}
): Promise<Response> => {
  const headers = {
    'Content-Type': 'application/json',
    ...options.headers,
  };

  return fetch(url, {
    ...options,
    headers,
    credentials: 'include',
  });
};

/**
 * Make a POST request
 */
export const post = async (url: string, data: any): Promise<Response> => {
  return authenticatedFetch(url, {
    method: 'POST',
    body: JSON.stringify(data),
  });
};

/**
 * Make a PUT request
 */
export const put = async (url: string, data: any): Promise<Response> => {
  return authenticatedFetch(url, {
    method: 'PUT',
    body: JSON.stringify(data),
  });
};

/**
 * Make a DELETE request
 */
export const del = async (url: string): Promise<Response> => {
  return authenticatedFetch(url, {
    method: 'DELETE',
  });
};

/**
 * Make a GET request
 */
export const get = async (url: string): Promise<Response> => {
  return authenticatedFetch(url, {
    method: 'GET',
  });
}; 